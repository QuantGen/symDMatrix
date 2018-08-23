#' A Matrix-Like Class to Represent a Symmetric Matrix Partitioned into
#' File-Backed Blocks.
#'
#' A `symDMatrix` is a symmetric matrix partitioned into file-backed blocks.
#' This approach allows for very large symmetric matrices, commonly found for
#' example when computing genetic relationship matrices on large cohorts. A
#' `symDMatrix` object behaves similarly to a regular `matrix` by implementing
#' key methods such as `[`, `dim`, and `dimnames`.
#'
#' The `symDMatrix` class is a [LinkedMatrix::RowLinkedMatrix] that nests
#' multiple [LinkedMatrix::ColumnLinkedMatrix] objects containing blocks of
#' type `ff_matrix`. Because the matrix is symmetric, only the diagonal and
#' upper-triangular blocks need to be stored, but for more efficient queries,
#' the lower-triangular blocks are virtual transposes of their diagonal
#' counterparts.
#'
#' @example man/examples/symDMatrix.R
#' @seealso [initialize()][initialize,symDMatrix-method()] to create a
#' `symDMatrix` object from scratch, or preferably, [as.symDMatrix()] to create
#' a `symDMatrix` object from other objects.
#' @aliases symDMatrix-class
#' @export symDMatrix
#' @exportClass symDMatrix
symDMatrix <- setClass("symDMatrix", contains = "RowLinkedMatrix")


#' Create a New symDMatrix Instance.
#'
#' This method is run when a [symDMatrix-class] object is created using
#' `symDMatrix(...)` or `new("symDMatrix", ...)`.
#'
#' Several structural checks are performed on the passed blocks: there must be
#' at least one block, the blocks must be of type
#' [LinkedMatrix::ColumnLinkedMatrix], and the number of blocks must be
#' consistent across the [LinkedMatrix::ColumnLinkedMatrix] objects. Each block
#' must inherit from `ff_matrix` and have the same number of rows or columns as
#' blocks in the same row or column, respectively. Non-final blocks have to be
#' square, unless if there is only a single block, in which case that block
#' also has to be square.
#'
#' @param .Object The [symDMatrix-class] instance to be initialized. This
#' argument is passed in by R and can be ignored, but still needs to be
#' documented.
#' @param ... [LinkedMatrix::ColumnLinkedMatrix] objects containing blocks that
#' inherit from `ff_matrix`.
#' @return A [symDMatrix-class] object.
#' @example man/examples/initialize.R
#' @seealso [as.symDMatrix()] to create a [symDMatrix-class] object from other
#' objects.
#' @export
setMethod("initialize", "symDMatrix", function(.Object, ...) {
    blocks <- list(...)
    nBlocks <- length(blocks)
    # Stop if there are no blocks
    if (nBlocks == 0L) {
        stop("there needs to be at least one block")
    }
    # Stop if blocks are not of type ColumnLinkedMatrix
    if (!all(sapply(blocks, class) == "ColumnLinkedMatrix")) {
        stop("blocks need to be of type ColumnLinkedMatrix")
    }
    # Stop if the number of nested blocks is inconsistent
    if (length(unique(sapply(blocks, LinkedMatrix::nNodes))) > 1L) {
        stop("number of nested blocks is inconsistent")
    }
    rowDims <- sapply(blocks, nrow)
    colDims <- rep(NA_integer_, nBlocks)
    for (rowIndex in seq_len(nBlocks)) {
        rowBlocks <- blocks[[rowIndex]]
        # Stop if nested blocks do not inherit from ff_matrix
        if (!all(sapply(rowBlocks, inherits, "ff_matrix"))) {
            stop("nested blocks need to inherit from ff_matrix")
        }
        # Stop if nested blocks do not have the same number of columns per column of blocks
        if (all(is.na(colDims))) {
            colDims[] <- sapply(rowBlocks, ncol)
        } else {
            if (!all(sapply(rowBlocks, ncol) == colDims)) {
                stop("all nested blocks need the same number of columns per column of blocks")
            }
        }
    }
    if (nBlocks > 1L) {
        # Stop if non-final block is not square
        if (any(rowDims[-length(rowDims)] != colDims[-length(colDims)])) {
            stop("non-final blocks need to be square")
        }
    } else {
        # Stop if first block is not square
        if (nrow(blocks[[1]][[1]]) != ncol(blocks[[1]][[1]])) {
            stop("the first block needs to be square")
        }
    }
    # Call RowLinkedMatrix constructor
    .Object <- callNextMethod(.Object, ...)
    return(.Object)
})


#' Load symDMatrix Objects from .RData Files.
#'
#' This function is similar to [base::load()], but it also initializes the
#' `ff_matrix` blocks in the [symDMatrix-class] object.
#'
#' @param file The name of an .RData file to be loaded.
#' @param readonly Set to TRUE to forbid writing to existing files.
#' @param envir The environment where to load the data.
#' @export
load.symDMatrix <- function(file, readonly = FALSE, envir = parent.frame()) {
    # Load data into new environment
    loadingEnv <- new.env()
    load(file = file, envir = loadingEnv)
    names <- ls(envir = loadingEnv)
    for (name in names) {
        object <- get(name, envir = loadingEnv)
        # Initialize blocks of symDMatrix objects
        if (class(object) == "symDMatrix") {
            nBlocks <- nBlocks(object)
            for (i in 1L:nBlocks) {
                for (j in 1L:nBlocks) {
                    object[[i]][[j]] <- initializeBlock(object[[i]][[j]], path = dirname(file), readonly = readonly)
                }
            }
        }
        # Assign object to envir
        assign(name, object, envir = envir)
    }
    message("Loaded objects: ", paste0(names, collapse = ", "))
}


initializeBlock <- function(x, ...) {
    UseMethod("initializeBlock")
}


# Absolute paths to ff files are not stored, so the ff objects have to be
# loaded from the same directory as the RData file.
initializeBlock.ff_matrix <- function(x, path, readonly = FALSE, ...) {
    # Store current working directory and set working directory to path
    cwd <- getwd()
    setwd(path)
    # Open ff object
    ff::open.ff(x, readonly = readonly)
    # Restore the working directory
    setwd(cwd)
    return(x)
}


initializeBlock.default <- function(x, ...) {
    return(x)
}


#' Return the Number of Column/Row Blocks of a symDMatrix Object.
#'
#' This function returns the number of row blocks the original matrix has been
#' partitioned into.
#'
#' @param x A [symDMatrix-class] object.
#' @return The number of column/row blocks of a [symDMatrix-class] object.
#' @example man/examples/nBlocks.R
#' @export
nBlocks <- function(x) {
    LinkedMatrix::nNodes(x)
}


#' Return the Block Size of a symDMatrix Object.
#'
#' This function returns the block size of a [symDMatrix-class] object.
#'
#' The last block of a column/row may be smaller than the other blocks. Its
#' size can be retrieved by setting `last` to `TRUE`.
#'
#' @param x A [symDMatrix-class] object.
#' @param last A boolean indicating whether to return the block size of the
#' last (`TRUE`) column/row block or any of the other blocks (`FALSE`,
#' default).
#' @return The block size of a [symDMatrix-class] object.
#' @example man/examples/blockSize.R
#' @export
blockSize <- function(x, last = FALSE) {
    row <- x[[1L]]
    if (last) {
        ncol(row[[LinkedMatrix::nNodes(row)]])
    } else {
        ncol(row[[1L]])
    }
}


#' Return the Block Structure of a symDMatrix Object.
#'
#' This function returns the block structure of a [symDMatrix-class] object and
#' can be useful when implementing custom indexing techniques.
#'
#' @param x A [symDMatrix-class] object.
#' @return A matrix with three columns: the block number, the start index and
#' the end index.
#' @export
blockIndex <- function(x) {
    nNodes <- LinkedMatrix::nNodes(x)
    index <- matrix(nrow = nNodes, ncol = 3L)
    colnames(index) <- c("block", "ini", "end")
    end <- 0L
    for (i in 1L:nNodes) {
        ini <- end + 1L
        end <- ini + nrow(x[[i]][[1L]]) - 1L
        index[i, ] <- c(i, ini, end)
    }
    return(index)
}


#' Coerce an Object to a symDMatrix Object.
#'
#' Coerce an object to a [symDMatrix-class] object.
#'
#' @param x A numeric matrix.
#' @param ... Additional arguments.
#' @return A [symDMatrix-class] object.
#' @seealso [as.symDMatrix.matrix()] to coerce a matrix or
#' [as.symDMatrix.character()] to coerce a vector of path names to a
#' [symDMatrix-class] object.
#' @export
as.symDMatrix <- function(x, ...) {
    UseMethod("as.symDMatrix")
}


#' Coerce a Matrix to a symDMatrix Object.
#'
#' This function creates a [symDMatrix-class] from a numeric matrix that is
#' assumed to be symmetric.
#'
#' The input matrix is broken into blocks and each block is stored as an
#' `ff_matrix` object. In addition, a metadata object called `symDMatrix.RData`
#' is created to allow for easy reloading of the [symDMatrix-class] object.
#'
#' @param x A symmetric numeric matrix.
#' @param blockSize The number of rows and columns of each block. If `NULL`, a
#' single block of the same dimensions as `x` will be created. Defaults to
#' 5000.
#' @param vmode The vmode used to store the data in the `ff` objects.
#' @param folderOut A name for a folder where to store the data of the
#' resulting [symDMatrix-class] object.
#' @param ... Additional arguments (currently unused).
#' @return A [symDMatrix-class] object.
#' @seealso [load.symDMatrix()] to reload the [symDMatrix-class] object.
#' @export
as.symDMatrix.matrix <- function(x, blockSize = 5000L, vmode = "double", folderOut = randomString(), ...) {
    n <- nrow(x)
    if (ncol(x) != n) {
        stop("x must be a square matrix")
    }
    if (file.exists(folderOut)) {
        stop(folderOut, " already exists")
    }
    dir.create(folderOut)
    # Determine number of blocks from block size
    nBlocks <- as.integer(ceiling(nrow(x) / blockSize))
    # Determine subjects of each block
    index <- matrix(data = integer(), nrow = nBlocks, ncol = 2L)
    index[1L, ] <- c(1L, min(n, blockSize))
    if (nBlocks > 1L) {
        for (i in 2L:nBlocks) {
            index[i, 1L] <- index[(i - 1L), 2L] + 1L
            index[i, 2L] <- min(index[i, 1L] + blockSize - 1L, n)
        }
    }
    # Create nested list
    args <- vector(mode = "list", length = nBlocks)
    for (rowIndex in 1L:nBlocks) {
        rowRanges <- seq(index[rowIndex, 1L], index[rowIndex, 2L])
        rowArgs <- vector(mode = "list", length = nBlocks)
        for (colIndex in 1L:nBlocks) {
            colRanges <- seq(index[colIndex, 1L], index[colIndex, 2L])
            blockName <- paste0("data_", padDigits(rowIndex, nBlocks), "_", padDigits(colIndex, nBlocks), ".bin")
            block <- ff::ff(dim = c(length(rowRanges), length(colRanges)), vmode = vmode, initdata = x[rowRanges, colRanges], filename = paste0(folderOut, "/", blockName), dimnames = list(rownames(x)[rowRanges], colnames(x)[colRanges]))
            # Change ff path to a relative one
            bit::physical(block)$filename <- blockName
            if (colIndex >= rowIndex) {
                rowArgs[[colIndex]] <- block
            } else {
                rowArgs[[colIndex]] <- ff::vt(args[[colIndex]][[rowIndex]])
            }
        }
        args[[rowIndex]] <- do.call(LinkedMatrix::ColumnLinkedMatrix, rowArgs)
    }
    # Create symDMatrix object from args
    symDMatrix <- do.call(symDMatrix, args)
    save(symDMatrix, file = paste0(folderOut, "/symDMatrix.RData"))
    return(symDMatrix)
}


#' Coerce a Character Vector to a symDMatrix Object.
#'
#' This function creates a [symDMatrix-class] object from a character vector of
#' path names to `RData` files, each containing exactly one `ff_matrix` object
#' that is used as a block, and is useful for distributed computing where each
#' block is processed on a different node.
#'
#' The `RData` files must be ordered by block: `G11, G12, G13, ..., G1q, G22,
#' G23, ..., G2q, ..., Gqq`. The matrix-like objects are initialized similarly
#' to [load.symDMatrix()].
#'
#' @param x A character vector with path names to `RData` files.
#' @param ... Additional arguments (currently unused).
#' @return A [symDMatrix-class] object.
#' @seealso [base::list.files()] to create a character vector of file paths
#' that match a certain pattern.
#' @export
as.symDMatrix.character <- function(x, ...) {
    nBlocks <- as.integer((-1L + sqrt(1L + 4L * 2L * length(x))) / 2L)
    args <- vector(mode = "list", length = nBlocks)
    counter <- 1L
    for (i in 1L:nBlocks) {
        rowArgs <- vector(mode = "list", length = nBlocks)
        for (j in 1L:nBlocks) {
            if (j >= i) {
                loadingEnv <- new.env()
                file <- x[[counter]]
                load(file = file, envir = loadingEnv)
                names <- ls(envir = loadingEnv)
                # Make sure that exactly one object inherits from ff
                inherits <- sapply(names, function(name) {
                    object <- get(name, envir = loadingEnv)
                    inherits(object, "ff_matrix")
                })
                if (sum(inherits) != 1L) {
                    stop("only one object per RData file can inherit from ff_matrix")
                }
                object <- get(names[which(inherits)], envir = loadingEnv)
                # Initialize the matrix-like object
                object <- initializeBlock(object, path = dirname(file))
                rowArgs[[j]] <- object
                counter <- counter + 1L
            } else {
                rowArgs[[j]] <- ff::vt(args[[j]][[i]])
            }
        }
        args[[i]] <- do.call(LinkedMatrix::ColumnLinkedMatrix, rowArgs)
    }
    # Create symDMatrix object from args
    symDMatrix <- do.call(symDMatrix, args)
    return(symDMatrix)
}


randomString <- function(n = 10L) {
    paste(sample(c(0L:9L, letters, LETTERS), size = n, replace = TRUE), collapse = "")
}


padDigits <- function(x, total) {
    formatC(x, width = as.integer(log10(total) + 1L), format = "d", flag = "0")
}
