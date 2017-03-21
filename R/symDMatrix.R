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
#' @export symDMatrix
#' @exportClass symDMatrix
symDMatrix <- setClass("symDMatrix", slots = c(data = "list", centers = "numeric", scales = "numeric"))


#' Creates a New symDMatrix Instance.
#'
#' This method is run when a [symDMatrix-class] object is created using
#' `symDMatrix(...)` or `new("symDMatrix", ...)`.
#'
#' @param .Object The [symDMatrix-class] instance to be initialized. This
#' argument is passed in by R and can be ignored, but still needs to be
#' documented.
#' @param data A nested list to fill the `@@data` slot of the
#' [symDMatrix-class] object.
#' @param centers A numeric vector to fill the `@@centers` slot of the
#' [symDMatrix-class] object.
#' @param scales A numeric vector to fill the `@@scales` slot of the
#' [symDMatrix-class] object.
setMethod("initialize", "symDMatrix", function(.Object, data, centers, scales) {
    nBlocks <- length(data)
    # Test that there is at least one block
    if (nBlocks == 0L) {
        stop("data needs to contain at least one block")
    }
    # Test that data has the right structure
    blocksPerRow <- sapply(data, length)
    if (!identical(blocksPerRow, seq(nBlocks, 1L))) {
        stop("data needs to be a nested list in the following structure: [[G11, G12, G13, ..., G1q], [G22, G23, ..., G2q], [...], [Gqq]]")
    }
    # Block-level tests
    rowDims <- rep(NA_integer_, nBlocks)
    colDims <- rep(NA_integer_, nBlocks)
    for (i in seq(0L, nBlocks - 1L)) {
        for (j in seq(1L, nBlocks - i)) {
            block <- data[[i + 1L]][[j]]
            # Test that all blocks are ff objects
            if (!inherits(block, "ff_matrix")) {
                stop("data: all blocks need to be ff_matrix objects")
            }
            # Test that all blocks per row have the same number of rows
            if (is.na(rowDims[i + 1L])) {
                rowDims[i + 1L] <- nrow(block)
            } else {
                if (nrow(block) != rowDims[i + 1L]) {
                    stop("data: all blocks per row need the same number of rows")
                }
            }
            # Test that all blocks per column have the same number of columns
            if (is.na(colDims[j + i])) {
                colDims[j + i] <- ncol(block)
            } else {
                if (ncol(block) != colDims[j + i]) {
                    stop("data: all blocks per column need the same number of columns")
                }
            }
        }
    }
    if (nBlocks > 1L) {
        # Test that non-final blocks are square
        if (any(rowDims[-length(rowDims)] != colDims[-length(colDims)])) {
            stop("data: non-final blocks need to be square")
        }
    } else {
        # Test that the first block is square
        if (nrow(data[[1]][[1]]) != ncol(data[[1]][[1]])) {
            stop("data: the first block needs to be square")
        }
    }
    .Object@data <- data
    .Object@centers <- centers
    .Object@scales <- scales
    return(.Object)
})


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

    nargs <- nargs()

    nX <- nrow(x)

    # Single Index: x[i]
    if (nargs == 2L && !missing(i) && missing(j)) {

        singleIndex <- TRUE

        n <- 1L
        p <- length(i)

        # Convert single index to multi index
        k <- i - 1L
        paired_i <- k %% nX
        paired_j <- as.integer(k / nX)
        paired_i <- paired_i + 1L
        paired_j <- paired_j + 1L

    # No index and multi Index: x[], or x[, ], and x[i, j], x[i, ], or x[, j]
    } else {

        singleIndex <- FALSE

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

        # Create all combinations of i and j
        paired_i <- rep(i, each = p)
        paired_j <- rep(j, times = n)

    }

    # Switch indices for combinations in which i is larger than j to redirect
    # queries to the lower triangle to the upper triangle
    switch <- paired_i > paired_j
    flip <- paired_i[switch]
    paired_i[switch] <- paired_j[switch]
    paired_j[switch] <- flip

    # Retrieve block size
    blockSize <- blockSize(x)

    # Create retrieval index
    row_blocks <- as.integer(ceiling(paired_i / blockSize))
    col_blocks <- as.integer(ceiling(paired_j / blockSize))
    local_i <- paired_i - (row_blocks - 1L) * blockSize
    local_j <- paired_j - (col_blocks - 1L) * blockSize

    # Initialize output matrix
    names <- names.symDMatrix(x)
    if (!is.null(names) && !singleIndex) {
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
        # Load genotypes of symDMatrix objects
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
#' @return The number of column/row blocks of a [symDMatrix-class] object.
#' @export
nBlocks <- function(x) {
    length(x@data[[1L]])
}


#' Returns the Block Size of a symDMatrix Object.
#'
#' The last block of a column/row may be smaller than the other blocks. Its
#' size can be retrieved by setting `last` to `TRUE`.
#'
#' @param x A [symDMatrix-class] object.
#' @param last A boolean indicating whether to return the block size of the
#' last (`TRUE`) column/row block or any of the other blocks (`FALSE`,
#' default).
#' @return The block size of a [symDMatrix-class] object.
#' @export
blockSize <- function(x, last = FALSE) {
    row <- x@data[[1L]]
    if (last) {
        ncol(row[[length(row)]])
    } else {
        ncol(row[[1L]])
    }
}


#' Returns the Block Structure of a symDMatrix Object.
#'
#' @param x A [symDMatrix-class] object.
#' @return A matrix.
#' @export
blockIndex <- function(x) {
    n <- length(x@data)
    index <- matrix(nrow = n, ncol = 3L)
    colnames(index) <- c("block", "ini", "end")
    end <- 0L
    for (i in 1L:n) {
        ini <- end + 1L
        end <- ini + nrow(x@data[[i]][[1L]]) - 1L
        index[i, ] <- c(i, ini, end)
    }
    return(index)
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


#' @export
as.symDMatrix.list <- function(x, centers = 0L, scales = 1L) {
    counter <- 1L
    nBlocks <- as.integer((-1L + sqrt(1L + 4L * 2L * length(x))) / 2L)
    dataList <- vector(mode = "list", length = nBlocks)
    for (i in 1L:nBlocks) {
        dataList[[i]] <- vector(mode = "list", length = nBlocks - i)
        for (j in i:nBlocks) {
            loadingEnv <- new.env()
            load(file = x[[counter]], envir = loadingEnv)
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


randomString <- function(n = 10L) {
    paste(sample(c(0L:9L, letters, LETTERS), size = n, replace = TRUE), collapse = "")
}
