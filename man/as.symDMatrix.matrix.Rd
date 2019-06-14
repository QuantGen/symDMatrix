\name{as.symDMatrix.matrix}
\alias{as.symDMatrix.matrix}
\title{Coerce a Matrix to a symDMatrix Object}
\description{
This function creates a \code{symDMatrix} from a numeric matrix that is assumed
to be symmetric.
}
\usage{
\method{as.symDMatrix}{matrix}(x, blockSize = 5000L, vmode = "double",
  folderOut = randomString(), ...)
}
\arguments{
  \item{x}{
    A symmetric numeric matrix.
  }
  \item{blockSize}{
    The number of rows and columns of each block. If \code{NULL}, a single
    block of the same dimensions as \code{x} will be created. Defaults to 5000.
  }
  \item{vmode}{
    The vmode used to store the data in the \code{ff} objects.
  }
  \item{folderOut}{
    A name for a folder where to store the data of the resulting
    \code{symDMatrix} object.
  }
  \item{...}{
    Additional arguments (currently unused).
  }
}
\details{
The input matrix is broken into blocks and each block is stored as an
\code{ff_matrix} object. In addition, a metadata object called
\code{symDMatrix.RData} is created to allow for easy reloading of the
\code{symDMatrix} object.
}
\value{
A \code{symDMatrix} object.
}
\seealso{
\code{\link{load.symDMatrix}} to reload the \code{symDMatrix} object.
}
