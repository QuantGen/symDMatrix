#' A Class to Represent a Symmetric Matrix Paritioned into Memory-Mapped Blocks.
#'
#' A `symDMatrix` is a symmetric matrix partitioned into memory-mapped blocks.
#' Because the matrix is symmetric, only the diagonal and upper-triangular
#' blocks are stored. Each block is an `ff` object.
#'
#' Internally, the blocks are organized as a nested list in the `@@data` slot,
#' each list element representing one row of the symmetric matrix. All blocks
#' except the ones in the last column/row are expected to have the same
#' dimensions.
#'
#' @exportClass symDMatrix
setClass("symDMatrix", slots = c(data = "list", centers = "numeric", scales = "numeric"))


#' A Constructor for Creating symDMatrix Objects.
#'
#' A [symDMatrix-class] is a symmetric matrix partitioned into memory-mapped
#' blocks. Because the matrix is symmetric, only the diagonal and
#' upper-triangular blocks are stored. Each block is an `ff` object.
#'
#' @param dataFiles A character vector with names of the `ff` files that
#' contain the data needed to create the object. The files must be ordered by
#' block, G11, G12, G13, ..., G1q, G22, G23, ..., Gqq.
#' @param centers (numeric) A vector storing the means used, if any, when
#' creating the symmetric matrix.
#' @param scales (numeric) A vector storing the standard deviations used, if
#' any, when creating the symmetric matrix.
#' @return A [symDMatrix-class] object.
#' @export
symDMatrix <- function(dataFiles, centers = 0L, scales = 1L) {
    counter <- 1L
    nBlocks <- as.integer((-1L + sqrt(1L + 4L * 2L * length(dataFiles))) / 2L)
    dataList <- vector(mode = "list", length = nBlocks)
    for (i in 1L:nBlocks) {
        dataList[[i]] <- vector(mode = "list", length = nBlocks - i)
        for (j in i:nBlocks) {
            loadingEnv <- new.env()
            load(file = dataFiles[[counter]], envir = loadingEnv)
            # TODO: Assumes only one object per data file
            objectName <- ls(envir = loadingEnv)[1L]
            object <- get(objectName, envir = loadingEnv)
            dataList[[i]][[j - i + 1L]] <- object
            counter <- counter + 1L
        }
    }
    G <- new("symDMatrix", data = dataList, centers = centers, scales = scales)
    return(G)
}


#' @export
is.matrix.symDMatrix <- function(x) {
    TRUE
}


#' @export
dim.symDMatrix <- function(x) {
    p <- sum(sapply(x@data[[1L]], ncol))
    c(p, p)
}


#' @export
length.symDMatrix <- function(x) {
    prod(dim(x))
}


names.symDMatrix <- function(x) {
    blockNames <- lapply(x@data[[1L]], function(block) {
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


#' @export
`[.symDMatrix` <- function(x, i, j, drop = TRUE) {

    nX <- nrow(x)
    pX <- ncol(x)

    if (missing(i)) {
        i <- 1L:nX
    } else if (typeof(i) == "logical") {
        i <- rep_len(i, nX)
        i <- which(i)
    } else if (typeof(i) == "character") {
        i <- match(i, rownames(x))
    } else if (typeof(i) == "double") {
        i <- as.integer(i)
    }
    if (missing(j)) {
        j <- 1L:pX
    } else if (typeof(j) == "logical") {
        j <- rep_len(j, pX)
        j <- which(j)
    } else if (typeof(j) == "character") {
        j <- match(j, colnames(x))
    } else if (typeof(j) == "double") {
        j <- as.integer(j)
    }

    n <- length(i)
    p <- length(j)

    # Retrieve block size
    blockSize <- blockSize(x)

    # Create all combinations of i and j and switch indices for combinations in
    # which i is larger than j to redirect queries to the lower triangle to the
    # upper triangle
    paired_i <- rep(i, each = p)
    paired_j <- rep(j, times = n)
    switch <- paired_i > paired_j
    flip <- paired_i[switch]
    paired_i[switch] <- paired_j[switch]
    paired_j[switch] <- flip

    # Create retrieval index
    row_blocks <- as.integer(ceiling(paired_i / blockSize))
    col_blocks <- as.integer(ceiling(paired_j / blockSize))
    local_i <- paired_i - (row_blocks - 1L) * blockSize
    local_j <- paired_j - (col_blocks - 1L) * blockSize

    # Initialize output matrix
    names <- names.symDMatrix(x)
    if (!is.null(names)) {
        dimnames <- list(names[i], names[j])
    } else {
        dimnames <- NULL
    }
    OUT <- matrix(data = double(), nrow = n, ncol = p, dimnames = dimnames)

    # Create output index
    out_i <- rep(1L:n, each = p)
    out_j <- rep(1L:p, times = n)

    # Retrieve elements by block
    for (row_block in unique(row_blocks)) {
        row_block_matches <- row_blocks == row_block
        for (col_block in unique(col_blocks[row_block_matches])) {
            cur_block <- row_block_matches & col_blocks == col_block
            block <- x@data[[row_block]][[col_block - row_block + 1L]]
            out_idx <- (out_j[cur_block] - 1L) * n + out_i[cur_block]
            local_idx <- (local_j[cur_block] - 1L) * nrow(block) + local_i[cur_block]
            OUT[out_idx] <- block[local_idx]
        }
    }

    if (drop == TRUE && (n == 1L || p == 1L)) {
        return(OUT[, ])
    } else {
        return(OUT)
    }

}


#' A Function to Load symDMatrix Objects into an R Session.
#'
#' Conceptually this function is similar to [base::load()]. However,
#' [load.symDMatrix()] also opens the connections to the `ff` files.
#'
#' @param file The name of an .RData file (created using [base::save()]).
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
            for (i in 1L:nBlocks) {
                for (j in i:nBlocks) {
                    ff::open.ff(object@data[[i]][[j - i + 1L]])
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


#' Determines the Number of Column/Row Blocks of a symDMatrix Object.
#'
#' @param x A [symDMatrix-class] object.
#' @export
nBlocks <- function(x) length(x@data[[1L]])


#' Returns the Column/Row Block Size of a symDMatrix Object.
#'
#' The last column/row block may be smaller.
#'
#' @param x A [symDMatrix-class] object.
#' @export
blockSize <- function(x) nrow(x@data[[1L]][[1L]])


#' Returns the Block Structure of a symDMatrix Object.
#'
#' @param x A [symDMatrix-class] object.
#' @export
blocks <- function(x) {
    n <- length(x@data)
    OUT <- matrix(nrow = n, ncol = 3L)
    colnames(OUT) <- c("block", "ini", "end")
    end <- 0L
    for (i in 1L:n) {
        ini <- end + 1L
        end <- ini + nrow(x@data[[i]][[1L]]) - 1L
        OUT[i, ] <- c(i, ini, end)
    }
    return(OUT)
}


#' Coerce an Object to a symDMatrix Object.
#'
#' @param x A numeric matrix.
#' @param ... Additional arguments.
#' @return A [symDMatrix-class] object.
#' @export
as.symDMatrix <- function(x, ...) {
    UseMethod("as.symDMatrix")
}


#' @rdname as.symDMatrix
#' @export
as.symDMatrix.matrix <- function(x, nBlocks = 3L, vmode = "double", folder = randomString(), saveRData = TRUE) {

    n <- nrow(x)

    if (ncol(x) != n) {
        stop("x must by a square matrix")
    }

    # Save current working directory before switching to destination path to
    # support relative paths in ff objects
    curDir <- getwd()
    dir.create(folder)
    setwd(folder)

    blockSize <- as.integer(ceiling(n / nBlocks))

    # Determe block size and subjects of each block
    index <- matrix(data = integer(), nrow = nBlocks, ncol = 3L)
    index[1L, ] <- c(1L, 1L, blockSize)
    if (nBlocks > 1L) {
        for (i in 2L:nBlocks) {
            index[i, 1L] <- i
            index[i, 2L] <- index[(i - 1L), 3L] + 1L
            index[i, 3L] <- min(index[i, 2L] + blockSize - 1L, n)
        }
    }

    dataList <- vector(mode = "list", length = nBlocks)
    ini <- 1L
    end <- 0L
    for (i in 1L:nBlocks) {
        rowIndex <- seq(index[i, 2L], index[i, 3L])
        dataList[[i]] <- vector(mode = "list", length = nBlocks - i)
        for (j in i:nBlocks) {
            colIndex <- seq(index[j, 2L], index[j, 3L])
            k <- j - i + 1L
            block <- ff::ff(dim = c(length(rowIndex), length(colIndex)), vmode = vmode,
                            initdata = x[rowIndex, colIndex],
                            filename = paste0("data_", i, "_", j, ".bin"))
            colnames(block) <- colnames(x)[colIndex]
            rownames(block) <- rownames(x)[rowIndex]
            bit::physical(block)$pattern <- "ff"
            bit::physical(block)$filename <- paste0("data_", i, "_", j, ".bin")
            dataList[[i]][[k]] <- block
        }
    }
    G <- new("symDMatrix", data = dataList, centers = 0L, scales = 0L)
    if (saveRData) {
        save(G, file = "G.RData")
    }

    # Restore working directory
    setwd(curDir)

    return(G)
}


randomString <- function(n = 10L) {
    paste(sample(c(0L:9L, letters, LETTERS), size = n, replace = TRUE), collapse = "")
}
