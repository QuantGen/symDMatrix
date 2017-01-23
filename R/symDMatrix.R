#' An S4 class that represents a memory-mapped symmetric matrix.
#'
#' @exportClass symDMatrix
setClass("symDMatrix", slots = c(data = "list", centers = "numeric", scales = "numeric"))


#' A constructor for creating \code{\linkS4class{symDMatrix}} objects.
#'
#' A \code{\linkS4class{symDMatrix}} is a memory-mapped symmetric matrix. The
#' data of the matrix is stored in \code{ff} files. The matrix is chunked into
#' blocks, the data for the diagonal and off-diagonal blocks (because the matrix
#' is symmetric only the upper-triangular blocks are stored) is stored in
#' \code{ff} files.
#'
#' @param dataFiles A character vector with names of the \code{ff} files that
#'   contain the data needed to create the object. The files must be order by
#'   block, G11, G12, G13, ..., G1q, G22, G23,...., Gqq.
#' @param centers (numeric) A vector storing the means used, if any, when
#'   creating G.
#' @param scales (numeric) A vector storing the standard deviations used, if
#'   any, when creating G.
#' @return A \code{\linkS4class{symDMatrix}} object.
#' @export
symDMatrix <- function(dataFiles, centers = 0, scales = 1) {
    counter <- 1
    nBlocks <- (-1 + sqrt(1 + 4 * 2 * length(dataFiles))) / 2
    dataList <- vector(mode = "list", length = nBlocks)
    for (i in 1:nBlocks) {
        dataList[[i]] <- vector(mode = "list", length = nBlocks - i)
        for (j in i:nBlocks) {
            loadingEnv <- new.env()
            load(file = dataFiles[[counter]], envir = loadingEnv)
            # TODO: Assumes only one object per data file
            objectName <- ls(envir = loadingEnv)[1]
            object <- get(objectName, envir = loadingEnv)
            dataList[[i]][[j - i + 1]] <- object
            counter <- counter + 1
        }
    }
    G <- new("symDMatrix", data = dataList, centers = centers, scales = scales)
    return(G)
}


#' @export
dim.symDMatrix <- function(x) {
    p <- sum(sapply(x@data[[1]], ncol))
    c(p, p)
}


#' @export
length.symDMatrix <- function(x) {
    prod(dim(x))
}


names.symDMatrix <- function(x) {
    blockNames <- lapply(x@data[[1]], function(block) {
        colnames(block)
    })
    isNULL <- sapply(blockNames, function(blockName) {
        is.null(blockName)
    })
    if (any(isNULL)) {
        NULL
    } else {
        unlist(blockNames)
    }
}


#' @export
dimnames.symDMatrix <- function(x) {
    names <- names.symDMatrix(x)
    if (is.null(names)) {
        NULL
    } else {
        list(names, names)
    }
}


diag.ff <- function(x) {
    n <- min(dim(x))
    out <- vector(mode = "double", length = n)
    for (i in 1:n) {
        out[i] <- x[i, i]
    }
    return(out)
}


diag.symDMatrix <- function(x) {
    out <- unlist(lapply(seq_len(nBlocks(x)), function(i) {
        block <- x@data[[i]][[1]]
        diag.ff(block)
    }))
    names(out) <- names.symDMatrix(x)
    return(out)
}


#' Extract the diagonal of a \code{\linkS4class{symDMatrix}} object.
#'
#' @inheritParams base::diag
#' @export
setMethod("diag", signature = "symDMatrix", definition = diag.symDMatrix)


#' Coerce a RAM numeric matrix (assumed to be symmetric) into a
#' \code{\linkS4class{symDMatrix}}
#'
#' @param x A numeric matrix.
#' @param nBlocks The number of column (also row) blocks to be used.
#' @param vmode The vmode used to store the data in the \code{ff} objects.
#' @param folder A name for a folder where to store the data of the resulting
#'   \code{\linkS4class{symDMatrix}}
#' @param saveRData If TRUE, the metadata (the \code{\linkS4class{symDMatrix}})
#'   is saved using the name G.RData.
#' @return A \code{\linkS4class{symDMatrix}} object.
#' @export
as.symDMatrix <- function(x, nBlocks = 3, vmode = "double", folder = randomString(), saveRData = TRUE) {
    n <- nrow(x)
    if (ncol(x) != n) {
        stop("x must by a square matrix")
    }

    tmpDir <- getwd()
    dir.create(folder)
    setwd(folder)

    blockSize <- ceiling(n/nBlocks)

    ## Determining block size and subjects in each block
    TMP <- matrix(nrow = nBlocks, ncol = 3)
    TMP[, 1] <- 1:nBlocks
    TMP[1, 2] <- 1
    TMP[1, 3] <- blockSize
    if (nBlocks > 1) {
        for (i in 2:nBlocks) {
            TMP[i, 2] <- TMP[(i - 1), 3] + 1
            TMP[i, 3] <- min(TMP[i, 2] + blockSize - 1, n)
        }
    }
    ini <- 1
    end <- 0
    DATA <- list()

    for (i in 1:nBlocks) {
        rowIndex <- eval(parse(text = paste0(TMP[i, 2], ":", TMP[i, 3])))
        DATA[[i]] <- list()
        for (j in i:nBlocks) {
            colIndex <- eval(parse(text = paste0(TMP[j, 2], ":", TMP[j, 3])))
            k <- j - i + 1
            DATA[[i]][[k]] <- ff::ff(dim = c(length(rowIndex), length(colIndex)), vmode = vmode,
                                     initdata = as.vector(x[rowIndex, colIndex]),
                                     filename = paste0("data_", i, "_", j, ".bin"))
            colnames(DATA[[i]][[k]]) <- colnames(x)[colIndex]
            rownames(DATA[[i]][[k]]) <- rownames(x)[rowIndex]
            bit::physical(DATA[[i]][[k]])$pattern <- "ff"
            bit::physical(DATA[[i]][[k]])$filename <- paste0("data_", i, "_", j, ".bin")
        }
    }
    G <- new("symDMatrix", data = DATA, centers = 0, scales = 0)
    if (saveRData) {
        save(G, file = "G.RData")
    }
    setwd(tmpDir)
    return(G)
}


subset.symDMatrix <- function(x, i, j, drop) {

    nX <- nrow(x)
    pX <- ncol(x)
    if (missing(i)) {
        i <- 1:nX
    }
    if (missing(j)) {
        j <- 1:pX
    }
    if (class(i) == "logical") {
        i <- rep_len(i, nX)
        i <- which(i)
    } else if (class(i) == "character") {
        i <- match(i, rownames(x))
    }
    if (class(j) == "logical") {
        j <- rep_len(j, pX)
        j <- which(j)
    } else if (class(j) == "character") {
        j <- match(j, colnames(x))
    }

    # Retrieve block size
    # TODO: do not assume that all blocks have the same size
    blockSize <- ncol(x@data[[1]][[1]])

    # Create all combinations of i and j and switch indices for combinations in
    # which i is larger than j to redirect queries to the lower triangle to the
    # upper triangle
    global.i <- rep(i, each = length(j))
    global.j <- rep(j, times = length(i))
    switch <- global.i > global.j
    flip <- global.i[switch]
    global.i[switch] <- global.j[switch]
    global.j[switch] <- flip

    # Create retrieval index
    row.chunks <- ceiling(global.i / blockSize)
    col.chunks <- ceiling(global.j / blockSize)
    local.i <- global.i - (row.chunks - 1) * blockSize
    local.j <- global.j - (col.chunks - 1) * blockSize

    # Initialize output matrix
    names <- names.symDMatrix(x)
    if (!is.null(names)) {
        dimnames <- list(names[i], names[j])
    } else {
        dimnames <- NULL
    }
    OUT <- matrix(data = double(), nrow = length(i), ncol = length(j), dimnames = dimnames)

    # Create output index
    out.i <- rep(1:length(i), each = length(j))
    out.j <- rep(1:length(j), times = length(i))

    # Retrieve elements by block
    for (row.chunk in unique(row.chunks)) {
        for (col.chunk in unique(col.chunks[which(row.chunks == row.chunk)])) {
            cur.chunk <- which(row.chunks == row.chunk & col.chunks == col.chunk)
            chunk.idx <- (local.j[cur.chunk] - 1) * nrow(x@data[[row.chunk]][[col.chunk - row.chunk + 1]]) + local.i[cur.chunk]
            out.idx <- (out.j[cur.chunk] - 1) * length(i) + out.i[cur.chunk]
            OUT[out.idx] <- x@data[[row.chunk]][[col.chunk - row.chunk + 1]][chunk.idx]
        }
    }

    if (drop == TRUE && (length(i) == 1 || length(j) == 1)) {
        return(OUT[, ])
    } else {
        return(OUT)
    }
}


#' Extract parts of a \code{\linkS4class{symDMatrix}}.
#'
#' @inheritParams base::`[`
#' @param j Column indices.
#' @export
setMethod("[", signature = "symDMatrix", definition = subset.symDMatrix)


#' A function to load \code{\linkS4class{symDMatrix}} objects into an R session.
#'
#' Conceptually this function is similar to \code{load()}. However,
#' \code{load.symDMatrix} also opens the connections to the \code{ff} files.
#'
#' @param file The name of an .RData file (created using \code{save()}).
#' @param envir The environment where to load the data.
#' @export
load.symDMatrix <- function(file, envir = parent.frame()) {
    # Load data into new environment
    loadingEnv <- new.env()
    load(file = file, envir = loadingEnv)
    names <- ls(envir = loadingEnv)
    for (name in names) {
        object <- get(name, envir = loadingEnv)
        # Load genotypes of BGData objects
        if (class(object) == "symDMatrix") {
            # Store current working directory and set working directory to
            # dirname of file
            cwd <- getwd()
            setwd(dirname(file))
            # Open ff objects
            nBlocks <- nBlocks(object)
            for (i in 1:nBlocks) {
                for (j in i:nBlocks) {
                    ff::open.ff(object@data[[i]][[j - i + 1]])
                }
            }
            # Restore the working directory
            setwd(cwd)
        }
        # Assign object to envir
        assign(name, object, envir = envir)
    }
    message("Loaded objects: ", paste0(names, collapse = ", "))
}


#' Determines the number of column/row blocks of a
#' \code{\linkS4class{symDMatrix}} object.
#'
#' @param x A \code{\linkS4class{symDMatrix}} object.
#' @export
nBlocks <- function(x) length(x@data[[1]])


#' Returns the column (also row) block size of a \code{\linkS4class{symDMatrix}}
#' object. Note, the last column/row block may be smaller.
#'
#' @param x A \code{\linkS4class{symDMatrix}} object.
#' @export
blockSize <- function(x) nrow(x@data[[1]][[1]])


#' Returns the block structure of a \code{\linkS4class{symDMatrix}} object.
#'
#' @param x A \code{\linkS4class{symDMatrix}} object.
#' @export
blocks <- function(x) {
    n <- length(x@data)
    OUT <- matrix(nrow = n, ncol = 3)
    OUT[, 1] <- 1:n
    colnames(OUT) <- c("block", "ini", "end")
    end <- 0
    for (i in 1:n) {
        ini <- end + 1
        end <- ini + nrow(x@data[[i]][[1]]) - 1
        OUT[i, 2] <- ini
        OUT[i, 3] <- end
        ini <- end + 1
    }
    return(OUT)
}


randomString <- function(n = 10) {
    paste(sample(c(0:9, letters, LETTERS), size = n, replace = TRUE), collapse = "")
}
