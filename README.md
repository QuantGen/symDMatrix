symDMatrix
==========

[![Travis-CI Build Status](https://travis-ci.org/QuantGen/symDMatrix.svg?branch=master)](https://travis-ci.org/QuantGen/symDMatrix)

symDMatrix is an R package that provides symmetric matrices assembled from memory-mapped blocks.

A symmetric matrix is partitioned into blocks as follows:

| G11 | G12 | G13 |
|:---:|:---:|:---:|
| G21 | G22 | G23 |
| G31 | G32 | G33 |

Because the matrix is assumed to be symmetric (i.e., Gij equals Gji), only the upper-triangular blocks are stored. Each block is stored as a flat file using an `ff` object.

The package defines the class and multiple methods that allow treating this memory-mapped matrix as a standard RAM matrix.

Internally, a `symDMatrix` object is an S4 class with the following slots:

* `@data` (list) each element of the list is an ff object
* `@centers` (numeric) column-means used in the computation of the matrix
* `@scales` (numeric) column-standard deviations used to scale the matrix


Tutorial
--------

### (0) Creating a symmetric matrix in RAM

Before we start, let's create a symmetric matrix in RAM.

```R
# Load genotypes from a mice data set
library(BGLR)
data(mice)

# Compute a symmetric genetic relationship matrix (G-matrix) in RAM
X <- mice.X
rownames(X) <- paste0("ID_", 1:nrow(X))
G <- tcrossprod(scale(X))
G <- G / mean(diag(G))
```

### (1) Converting a RAM matrix into a symDMatrix

In practice, if we can hold a matrix in RAM, there is not much of a point to convert it to a `symDMatrix` object; however, this will help us to get started.

```R
library(symDMatrix)

G2 <- as.symDMatrix(G, folder = "mice", nBlocks = 5, vmode = "double") # can use `single` for lighter files
```

### (2) Exploring operators

Now that we have a `symDMatrix` object, let's illustrate some operators.

```R
# Basic operators applied to a matrix in RAM and to a symDMatrix

# Dimension operators
dim(G) == dim(G2)
nrow(G) == nrow(G2)
ncol(G) == ncol(G2)
all.equal(diag(G), diag(G2))

# Names operators
all(dimnames(G)[[1]] == dimnames(G2)[[1]])
all(dimnames(G)[[2]] == dimnames(G2)[[2]])
all(rownames(G) == rownames(G2))
all(colnames(G) == rownames(G2))

# Block operators
nBlocks(G2)
blockSize(G2)
blocks(G2)

# Indexing (can use booleans, integers or labels)
G2[1:2, 1:2]
G2[c("ID_1", "ID_2"), c("ID_1", "ID_2")]
tmp <- c(TRUE, TRUE, rep(FALSE, nrow(G2) - 2))
G2[tmp, tmp]
head(G2[tmp, ])

# Exhaustive check of indexing
for (i in 1:100) {
    n1 <- sample(1:50, size = 1)
    n2 <- sample(1:50, size = 1)
    i1 <- sample(1:nrow(X), size = n1)
    i2 <- sample(1:nrow(X), size = n2)
    TMP1 <- G[i1, i2, drop = FALSE]
    TMP2 <- G2[i1, i2, drop = FALSE]
    stopifnot(all.equal(TMP1, TMP2))
}

```

### (3) Creating a symDMatrix from genotypes

The function `getG.symDMatrix` of the [BGData](https://github.com/QuantGen/BGData) package computes G=XX' (with options for centering and scaling) without ever loading G in RAM. It creates the `symDMatrix` object directly. In this example, X is a matrix in RAM. For large genotype data sets, X could be a memory-mapped matrix, `ff` object, or part of a `BGData` object.

```R
library(BGData)

G3 <- getG.symDMatrix(X, scaleCol = TRUE, centerCol = TRUE, nBlocks = 5, folder = "tmp", vmode = "double")
class(G3)
all.equal(diag(G), diag(G3))

for(i in 1:10){
    n1 <- sample(1:25, size = 1)
    i1 <- sample(1:25, size = n1)
    for(j in 1:10){
        n2 <- sample(1:nrow(G), size = 1)
        i2 <- sample(1:nrow(G), size = n2)
        tmp1 <- G[i1, i2]
        tmp2 <- G3[i1, i2]
        stopifnot(all.equal(tmp1, tmp2))
    }
}
```

### (4) Creating a symDMatrix from `ff` files containing the blocks

The function `symDMatrix` allows creating a `symDMatrix` object from a list of `ff` files. The list is assumed to provide, in order, files for G11, G12, ..., G1q, G22, G23, ..., G2q, ..., Gqq. This approach will be useful for very large G matrices. If n is large it may make sense to compute the blocks of the `symDMatrix` object in parallel jobs (e.g., in an HPC). The function `getGij` is similar to `getG.symDMatrix` (see [BGData](https://github.com/QuantGen/BGData) package) but accepts arguments `i1` and `i2` which define a block of G (i.e., rows of X).

```R
library(BGData)
library(BGLR)
data(wheat)
X <- wheat.X
rownames(X) <- 1:nrow(X)

centers <- colMeans(X)
scales <- apply(X = X, MARGIN = 2, FUN = sd)

nBlocks <- 8

dir.create("GMatrix")
setwd("GMatrix")

stepSize <- ceiling(nrow(X) / nBlocks)

G1 <- tcrossprod(scale(X, center = centers, scale = scales))
G1 <- G1 / mean(diag(G1))

# This loop may be executed in parallel
for (i in 1:nBlocks) {
    i_ini <- (i - 1) * stepSize + 1
    if (i_ini <= nrow(X)) {
        i_end <- min(nrow(X), i_ini + stepSize - 1)
        i1<-i_ini:i_end
        for (j in i:nBlocks) {
            j_ini <- (j - 1) * stepSize + 1
            j_end <- min(nrow(X), j_ini + stepSize - 1)
            i2 <- j_ini:j_end
            if (j_ini <= nrow(X)) {
                Gij <- getG(x = X, i = i1, i2 = i2,
                            saveName = paste0("G_", i, "_", j, ".bin"),
                            saveType = "ff", centers = centers,
                            scales = scales, scaleG = TRUE, verbose = FALSE,
                            scaleCol = TRUE, nChunks = 4, nChunks2 = 4,
                            mc.cores = 4, saveG = TRUE)
                print(paste(i_ini, i_end, " ; ",  j_ini, j_end))
            }
        }
    }
}

# Note: file list needs to be orderd G11, G12, ..., Gqq
G2 <- symDMatrix(dataFiles = list.files(pattern = "*.ff"), names = rownames(X))
```


Installation
------------

To get the current released version from CRAN:

```r
install.packages("symDMatrix")
```

To get the current development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("QuantGen/symDMatrix")
```


Example Dataset
---------------

The example dataset in the `inst/extdata` folder is the G matrix of the dummy dataset that comes with the [BEDMatrix](https://cran.r-project.org/package=BEDMatrix) package. It has been generated as follows:

```R
library(BGData)

X <- BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix"))

G <- getG.symDMatrix(X, nBlocks = 3, folder = "inst/extdata")
```

To load the dataset:

```R
load.symDMatrix(system.file("extdata", "G.RData", package = "symDMatrix")) # loads G
```
