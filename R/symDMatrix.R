#' An S4 class that represents a memory-mapped symmetric matrix.
#'
#' @exportClass symDMatrix
setClass("symDMatrix", slots = c(names = "character", centers = "numeric", scales = "numeric", data = "list"))

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
#' @param names (character) The rownames of the matrix.
#' @return A \code{\linkS4class{symDMatrix}} object.
#' @export
symDMatrix <- function(dataFiles, centers = 0, scales = 1, names = character()) {
    if (is.list(dataFiles)) {
        dataFiles <- unlist(dataFiles)
    }
    counter <- 1
    dataList <- list()
    n_chunks <- (-1 + sqrt(1 + 4 * 2 * length(dataFiles)))/2
    for (i in 1:n_chunks) {
        dataList[[i]] <- list()
        for (j in i:n_chunks) {
            oldList <- ls()
            load(dataFiles[counter])
            newList <- ls()
            objectName <- newList[which(!newList %in% c("oldList", oldList))]
            dataList[[i]][[j - i + 1]] <- get(objectName)
            counter <- counter + 1
            rm(list = objectName)
        }
    }
    G <- new("symDMatrix", names = names, centers = centers, scales = scales, data = dataList)
    return(G)
}

#' @export
dim.symDMatrix <- function(x) {
    rep(length(x@names), 2)
}

rownames.symDMatrix <- function(x) x@names

colnames.symDMatrix <- function(x) x@names

#' @export
dimnames.symDMatrix <- function(x) list(rownames.symDMatrix(x), colnames.symDMatrix(x))

#' @export
`dimnames<-.symDMatrix` <- function(x, value) {
    x@names <- value[[1]]
    return(x)
}

diag.ff <- function(x) {
    if (class(x)[1] != "ff_matrix") {
        stop("x must be an ff_matrix object")
    }
    n <- min(dim(x))
    out <- rep(NA, n)
    for (i in 1:n) {
        out[i] <- x[i, i]
    }
    return(out)
}

diag.symDMatrix <- function(x) {
    n <- min(dim(x))
    out <- rep(NA, n)

    nChunks <- nChunks(x)
    end <- 0
    for (i in 1:nChunks) {
        tmp <- diag.ff(x@data[[i]][[1]])
        ini <- end + 1
        end <- ini + length(tmp) - 1
        out[ini:end] <- tmp
    }
    names(out) <- x@names
    return(out)
}

#' Extract the diagonal of a \code{\linkS4class{symDMatrix}}
#'
#' @inheritParams base::diag
#' @export
setMethod("diag", signature = "symDMatrix", definition = diag.symDMatrix)


#' Coerce a RAM numeric matrix (assumed to be symmetric) into a
#' \code{\linkS4class{symDMatrix}}
#'
#' @param x A numeric matrix.
#' @param nChunks The number of column (also row) blocks to be used.
#' @param vmode The vmode used to store the data in the \code{ff} objects.
#' @param folder A name for a folder where to store the data of the resulting
#'   \code{\linkS4class{symDMatrix}}
#' @param saveRData If TRUE, the metadata (the \code{\linkS4class{symDMatrix}})
#'   is saved using the name G.RData.
#' @return A \code{\linkS4class{symDMatrix}} object.
#' @export
as.symDMatrix <- function(x, nChunks = 3, vmode = "double", folder = randomString(), saveRData = TRUE) {
    n <- nrow(x)
    if (ncol(x) != n) {
        stop("x must by a square matrix")
    }

    tmpDir <- getwd()
    dir.create(folder)
    setwd(folder)

    chunkSize <- ceiling(n/nChunks)

    ## Determining chunk size and subjects in each chunk
    TMP <- matrix(nrow = nChunks, ncol = 3)
    TMP[, 1] <- 1:nChunks
    TMP[1, 2] <- 1
    TMP[1, 3] <- chunkSize
    if (nChunks > 1) {
        for (i in 2:nChunks) {
            TMP[i, 2] <- TMP[(i - 1), 3] + 1
            TMP[i, 3] <- min(TMP[i, 2] + chunkSize - 1, n)
        }
    }
    ini <- 1
    end <- 0
    DATA <- list()

    for (i in 1:nChunks) {
        rowIndex <- eval(parse(text = paste0(TMP[i, 2], ":", TMP[i, 3])))
        DATA[[i]] <- list()
        for (j in i:nChunks) {
            colIndex <- eval(parse(text = paste0(TMP[j, 2], ":", TMP[j, 3])))
            k <- j - i + 1
            DATA[[i]][[k]] <- ff(dim = c(length(rowIndex), length(colIndex)), vmode = vmode,
                                 initdata = as.vector(x[rowIndex, colIndex]), filename = paste0("data_",
                                 i, "_", j, ".bin"))
            colnames(DATA[[i]][[k]]) <- colnames(x)[colIndex]
            rownames(DATA[[i]][[k]]) <- rownames(x)[rowIndex]
            physical(DATA[[i]][[k]])$pattern <- "ff"
            physical(DATA[[i]][[k]])$filename <- paste0("data_", i, "_", j, ".bin")
        }
    }
    tmp <- rownames(x)
    if (is.null(tmp)) {
        tmp <- paste0("id_", 1:nrow(x))
    }
    G <- new("symDMatrix", names = tmp, data = DATA, centers = 0, scales = 0)
    if (saveRData) {
        save(G, file = "G.RData")
    }
    setwd(tmpDir)
    return(G)
}

subset.symDMatrix <- function(x, i, j, drop) {

    if (missing(i)) {
        i <- 1:nrow(x)
    }
    if (missing(j)) {
        j <- 1:ncol(x)
    }
    if (class(i) == "logical") {
        i <- which(i)
    } else if (class(i) == "character") {
        i <- sapply(i, function(name) {
            which(rownames(x) == name)
        }, USE.NAMES = FALSE)
    }
    if (class(j) == "logical") {
        j <- which(j)
    } else if (class(j) == "character") {
        j <- sapply(j, function(name) {
            which(colnames(x) == name)
        }, USE.NAMES = FALSE)
    }

    nChunks <- nChunks(x)
    chunkSize <- ncol(x@data[[1]][[1]])
    i0 <- i
    j0 <- j
    i <- rep(i0, each = length(j0))
    j <- rep(j0, times = length(i0))

    # switching indexes ( i must be smaller than j)
    tmp <- i > j
    tmp2 <- i[tmp]
    i[tmp] <- j[tmp]
    j[tmp] <- tmp2

    out.i <- rep(1:length(i0), each = length(j0))
    out.j <- rep(1:length(j0), length(i0))

    row.chunk <- ceiling(i/chunkSize)
    col.chunk <- ceiling(j/chunkSize)
    local.i <- i - (row.chunk - 1) * chunkSize
    local.j <- j - (col.chunk - 1) * chunkSize

    OUT <- matrix(nrow = length(i0), ncol = length(j0), NA)
    rownames(OUT) <- x@names[i0]
    colnames(OUT) <- x@names[j0]

    for (i in unique(row.chunk)) {
        tmp <- which(row.chunk == i)
        for (j in unique(col.chunk[tmp])) {
            k <- which(row.chunk == i & col.chunk == j)
            tmp.row.in <- local.i[k]
            tmp.col.in <- local.j[k]
            tmp.in <- (tmp.col.in - 1) * nrow(x@data[[i]][[j - i + 1]]) + tmp.row.in
            tmp.row.out <- out.i[k]
            tmp.col.out <- out.j[k]
            tmp.out <- (tmp.col.out - 1) * nrow(OUT) + tmp.row.out
            OUT[tmp.out] <- x@data[[i]][[j - i + 1]][tmp.in]
        }
    }

    if (drop == TRUE && (length(i0) == 1 || length(j0) == 1)) {
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
#' @param verbose Whether to print additional information.
#' @export
load.symDMatrix <- function(file, envir = parent.frame(), verbose = TRUE) {
    # determining the object name
    lsOLD <- ls()
    load(file = file)
    lsNEW <- ls()
    objectName <- lsNEW[(!lsNEW %in% lsOLD) & (lsNEW != "lsOLD")]

    # determining path and filename
    path <- dirname(file)
    fname <- basename(file)

    # stores current working directiory and sets working directory to path
    cwd <- getwd()
    setwd(path)

    # determining object class
    objectClass <- class(eval(parse(text = objectName)))

    if (verbose) {
        cat(" Meta data (", fname, ") and its data were stored at folder ", path, ".\n", sep = "")
        cat(" Object Name: ", objectName, "\n", sep = "")
        cat(" Object Class: ", objectClass, "\n", sep = "")
    }
    if (!(objectClass %in% c("BGData", "rmmMatrix", "cmmMatrix", "symDMatrix"))) {
        stop(" Object class must be either BGData, cmmMatrix, rmmMatrix or symDMatrix")
    }

    # Determining number of chunks
    nChunks <- nChunks(eval(parse(text = objectName)))

    # opening files
    for (i in 1:nChunks) {
        for (j in i:nChunks) {
            if (verbose) {
                cat(" Opening flat file ", i, "\n")
            }
            open(eval(parse(text = paste0(objectName, "@data[[", i, "]][[", j - i + 1, "]]"))))
        }
    }
    # sending the object to envir
    assign(objectName, get(objectName), envir = envir)

    # restoring the working directory
    setwd(cwd)
    if (verbose) {
        cat(" Original directory (", getwd(), ") restored \n", sep = "")
    }
}

#' Determines the number of column/row chunks of a
#' \code{\linkS4class{symDMatrix}} object.
#'
#' @param x A \code{\linkS4class{symDMatrix}} object.
#' @export
nChunks <- function(x) length(x@data[[1]])

#' Returns the column (also row) chunk size of a \code{\linkS4class{symDMatrix}}
#' object. Note, the last column/row block may be smaller.
#'
#' @param x A \code{\linkS4class{symDMatrix}} object.
#' @export
chunkSize <- function(x) nrow(x@data[[1]][[1]])

#' Returns the chunk structure of a \code{\linkS4class{symDMatrix}} object.
#'
#' @param x A \code{\linkS4class{symDMatrix}} object.
#' @export
chunks <- function(x) {
    if (class(x) != "symDMatrix") {
        stop(" the input must be a symDMatrix object.")
    }

    n <- length(x@data)
    OUT <- matrix(nrow = n, ncol = 3)
    OUT[, 1] <- 1:n
    colnames(OUT) <- c("chunk", "ini", "end")
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

randomString <- function(n = 10) paste(sample(c(0:9, letters, LETTERS), size = n, replace = TRUE), collapse = "")
