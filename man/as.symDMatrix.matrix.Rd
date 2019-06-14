\name{as.symDMatrix.matrix}
\alias{as.symDMatrix.matrix}
\title{Coerce a Matrix to a symDMatrix Object.}
\usage{
\method{as.symDMatrix}{matrix}(x, blockSize = 5000L, vmode = "double",
  folderOut = randomString(), ...)
}
\arguments{
\item{x}{A symmetric numeric matrix.}

\item{blockSize}{The number of rows and columns of each block. If \code{NULL}, a
single block of the same dimensions as \code{x} will be created. Defaults to
5000.}

\item{vmode}{The vmode used to store the data in the \code{ff} objects.}

\item{folderOut}{A name for a folder where to store the data of the
resulting \linkS4class{symDMatrix} object.}

\item{...}{Additional arguments (currently unused).}
}
\value{
A \linkS4class{symDMatrix} object.
}
\description{
This function creates a \linkS4class{symDMatrix} from a numeric matrix that is
assumed to be symmetric.
}
\details{
The input matrix is broken into blocks and each block is stored as an
\code{ff_matrix} object. In addition, a metadata object called \code{symDMatrix.RData}
is created to allow for easy reloading of the \linkS4class{symDMatrix} object.
}
\seealso{
\code{\link[=load.symDMatrix]{load.symDMatrix()}} to reload the \linkS4class{symDMatrix} object.
}
