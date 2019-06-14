symDMatrix <- setClass("symDMatrix", contains = "RowLinkedMatrix")

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

nBlocks <- function(x) {
    LinkedMatrix::nNodes(x)
}

blockSize <- function(x, last = FALSE) {
    row <- x[[1L]]
    if (last) {
        ncol(row[[LinkedMatrix::nNodes(row)]])
    } else {
        ncol(row[[1L]])
    }
}

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

as.symDMatrix <- function(x, ...) {
    UseMethod("as.symDMatrix")
}

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
