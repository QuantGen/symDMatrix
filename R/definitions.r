setOldClass("ff_matrix")

`colnames<-.symDMatrix` <- function(x, value) {
    x@names <- values
    return(x)
}

`rownames<-.symDMatrix` <- function(x, value) {
    x@names <- values
    return(x)
}

`rownames<-.symDMatrix` <- function(x, value) {
    x@names <- values
    return(x)
}

setClass("symDMatrix", slots = c(names = "character", centers = "numeric", scales = "numeric", 
    data = "list"))

# An interface for creating symDMatrix objects
symDMatrix <- function(dataFiles, centers = 0, scales = 1, names = character()) {
    if (is.list(fileList)) {
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

nChunks <- function(x) length(x@data[[1]])
chunkSize <- function(x) nrow(x@data[[1]][[1]])

setMethod("rownames", signature = "symDMatrix", definition = function(x) x@names)
setMethod("colnames", signature = "symDMatrix", definition = function(x) x@names)
setMethod("dimnames", signature = "symDMatrix", definition = function(x) list(rownames(x), 
    colnames(x)))
nrow.symDMatrix <- function(x) length(x@names)
setMethod("nrow", signature = "symDMatrix", definition = nrow.symDMatrix)
setMethod("ncol", signature = "symDMatrix", definition = nrow.symDMatrix)
setMethod("dim", signature = "symDMatrix", definition = function(x) rep(nrow(x), 
    2))

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
# setMethod('diag',signature='ff_matrix',definition=diag.ff)

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
setMethod("diag", signature = "symDMatrix", definition = diag.symDMatrix)

as.symDMatrix <- function(x, nChunks = 3, vmode = "double", folder = randomString(), 
    saveRData = TRUE) {
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
    return(OUT)
}
setMethod("[", signature = "symDMatrix", definition = subset.symDMatrix)


getG.symDMatrix <- function(X, nChunks = 5, chunkSize = NULL, centers = NULL, scales = NULL, 
    centerCol = T, scaleCol = T, nChunks2 = 1, folder = randomString(5), vmode = "double", 
    verbose = TRUE, saveRData = TRUE, mc.cores = 1, scaleG = T) {
    
    timeIn <- proc.time()[3]
    n <- nrow(X)
    p <- ncol(X)
    if (is.null(chunkSize)) {
        chunkSize <- ceiling(n/nChunks)
    }
    
    if ((centerCol | scaleCol) & (is.null(centers) | is.null(scales))) {
        if (is.null(centers) & is.null(scales)) {
            centers <- rep(NA, p)
            scales <- rep(NA, p)
            for (i in 1:p) {
                xi <- X[, i]
                scales[i] <- sd(xi, na.rm = TRUE) * sqrt((n - 1)/n)
                centers[i] <- mean(xi, na.rm = TRUE)
            }
        }
        if ((!is.null(centers)) & (is.null(scales))) {
            scales <- rep(NA, p)
            for (i in 1:p) {
                xi <- X[, i]
                scales[i] <- sd(xi, na.rm = TRUE) * sqrt((n - 1)/n)
            }
        }
        if ((is.null(centers)) & (!is.null(scales))) {
            centers <- rep(NA, p)
            for (i in 1:p) {
                xi <- X[, i]
                centers[i] <- mean(xi, na.rm = TRUE)
            }
        }
    }
    
    if (!centerCol) 
        centers <- rep(0, p)
    if (!scaleCol) 
        scales <- rep(1, p)
    
    chunkID <- ceiling(1:n/chunkSize)
    nChunks <- max(chunkID)
    nFiles <- nChunks * (nChunks + 1)/2
    DATA <- list()
    counter <- 1
    
    tmpDir <- getwd()
    dir.create(folder)
    setwd(folder)
    
    for (i in 1:nChunks) {
        DATA[[i]] <- list()
        rowIndex_i <- which(chunkID == i)
        Xi <- X[rowIndex_i, ]
        
        # centering/scaling
        for (k in 1:p) {
            xik <- Xi[, k]
            xik <- (xik - centers[k])/scales[k]
            xik[is.na(xik)] <- 0
            Xi[, k] <- xik
        }
        
        for (j in i:nChunks) {
            rowIndex_j <- which(chunkID == j)
            Xj <- X[rowIndex_j, ]
            
            # centering/scaling
            for (k in 1:p) {
                xjk <- Xj[, k]
                xjk <- (xjk - centers[k])/scales[k]
                xjk[is.na(xjk)] <- 0
                Xj[, k] <- xjk
            }
            
            Gij <- tcrossprod.parallel(x = Xi, y = Xj, mc.cores = mc.cores, nChunks = nChunks2)
            
            DATA[[i]][[j - i + 1]] <- ff(dim = dim(Gij), vmode = vmode, initdata = as.vector(Gij), 
                filename = paste0("data_", i, "_", j, ".bin"))
            colnames(DATA[[i]][[j - i + 1]]) <- colnames(X)[rowIndex_j]
            rownames(DATA[[i]][[j - i + 1]]) <- rownames(X)[rowIndex_i]
            counter <- counter + 1
            physical(DATA[[i]][[j - i + 1]])$pattern <- "ff"
            physical(DATA[[i]][[j - i + 1]])$filename <- paste0("data_", i, "_", 
                j, ".bin")
            
            if (verbose) {
                cat(" Done with pair ", i, "-", j, " (", round(100 * counter/(nChunks * 
                  (nChunks + 1)/2)), "% ", round(proc.time()[3] - timeIn, 3), " seconds).\n", 
                  sep = "")
            }
        }
    }
    if (is.null(rownames(X))) 
        rownames(X) <- 1:n
    names(centers) <- colnames(X)
    names(scales) <- colnames(X)
    G <- new("symDMatrix", names = rownames(X), data = DATA, centers = centers, scales = scales)
    if (scaleG) {
        K <- mean(diag(G))
        for (i in 1:length(G@data)) {
            for (j in 1:length(G@data[[i]])) {
                G@data[[i]][[j]][] <- G@data[[i]][[j]][]/K
            }
        }
    }
    if (saveRData) {
        save(G, file = "G.RData")
    }
    setwd(tmpDir)
    return(G)
}



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
        cat(" Meta data (", fname, ") and its data were stored at folder ", path, 
            ".\n", sep = "")
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
            open(eval(parse(text = paste0(objectName, "@data[[", i, "]][[", j - i + 
                1, "]]"))))
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

randomString <- function(n = 10) paste(sample(c(0:9, letters, LETTERS), size = n, 
    replace = TRUE), collapse = "")
 
