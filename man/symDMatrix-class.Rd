\docType{class}
\name{symDMatrix-class}
\alias{symDMatrix-class}
\title{A Matrix-Like Class to Represent a Symmetric Matrix Partitioned into
    File-Backed Blocks}
\description{
    A \code{symDMatrix} is a symmetric matrix partitioned into file-backed
    blocks.  This approach allows for very large symmetric matrices, commonly
    found for example when computing genetic relationship matrices on large
    cohorts. A \code{symDMatrix} object behaves similarly to a regular
    \code{matrix} by implementing key methods such as \code{[}, \code{dim}, and
    \code{dimnames}.
}
\details{
    The \code{symDMatrix} class is a \code{RowLinkedMatrix} that nests multiple
    \code{ColumnLinkedMatrix} objects containing blocks of type
    \code{ff_matrix}.  Because the matrix is symmetric, only the diagonal and
    upper-triangular blocks need to be stored, but for more efficient queries,
    the lower-triangular blocks are virtual transposes of their diagonal
    counterparts.
}
\seealso{
    \code{\link{symDMatrix}} to create a \code{symDMatrix} object from scratch,
    or preferably, \code{\link{as.symDMatrix}} to create a \code{symDMatrix}
    object from other objects.
}
\examples{
# Get the path to the example symmetric matrix
path <- system.file("extdata", "G.RData", package = "symDMatrix")

# Load the example symDMatrix object (G)
load.symDMatrix(path, readonly = TRUE)

# Get the dimensions
dim(G)

# Get the row names
rownames(G)

# Get the column names
colnames(G)

# Extract the diagonal
diag(G)

# Extract rows and columns
G[1, ]
G[1:3, ]
G["per0_per0", ]
G[c("per0_per0", "per1_per1", "per2_per2"), ]
}
