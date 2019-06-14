\docType{methods}
\name{initialize,symDMatrix-method}
\alias{initialize,symDMatrix-method}
\title{Create a New symDMatrix Instance}
\description{
This method is run when a \linkS4class{symDMatrix} object is created using
\code{symDMatrix(...)} or \code{new("symDMatrix", ...)}.
}
\usage{
\S4method{initialize}{symDMatrix}(.Object, ...)
}
\arguments{
\item{.Object}{The \linkS4class{symDMatrix} instance to be initialized. This
argument is passed in by R and can be ignored, but still needs to be
documented.}

\item{...}{\link[LinkedMatrix:ColumnLinkedMatrix]{LinkedMatrix::ColumnLinkedMatrix} objects containing blocks that
inherit from \code{ff_matrix}.}
}
\details{
Several structural checks are performed on the passed blocks: there must be
at least one block, the blocks must be of type
\link[LinkedMatrix:ColumnLinkedMatrix]{LinkedMatrix::ColumnLinkedMatrix}, and the number of blocks must be
consistent across the \link[LinkedMatrix:ColumnLinkedMatrix]{LinkedMatrix::ColumnLinkedMatrix} objects. Each block
must inherit from \code{ff_matrix} and have the same number of rows or columns as
blocks in the same row or column, respectively. Non-final blocks have to be
square, unless if there is only a single block, in which case that block
also has to be square.
}
\value{
A \linkS4class{symDMatrix} object.
}
\seealso{
\code{\link[=as.symDMatrix]{as.symDMatrix()}} to create a \linkS4class{symDMatrix} object from other
objects.
}
\examples{
# Generate a symmetric matrix
X <- cov(matrix(data = rnorm(25), nrow = 5, ncol = 5))

# Break this matrix into blocks X11, X12, X22
# X21 can be stored as a virtual transpose of X12
X11 <- ff::as.ff(X[1:3, 1:3])
X12 <- ff::as.ff(X[1:3, 4:5])
X22 <- ff::as.ff(X[4:5, 4:5])
X21 <- ff::vt(X12)

# Create a symDMatrix from blocks
S <- symDMatrix(
    LinkedMatrix::ColumnLinkedMatrix(X11, X12),
    LinkedMatrix::ColumnLinkedMatrix(X21, X22)
)
nBlocks(S)
blockSize(S)
blockSize(S, last = TRUE)
}
