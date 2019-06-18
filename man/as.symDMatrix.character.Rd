\name{as.symDMatrix.character}
\alias{as.symDMatrix.character}
\title{Coerce a Character Vector to a symDMatrix Object}
\description{
    This function creates a \code{symDMatrix} object from a character vector of
    path names to \code{RData} files, each containing exactly one
    \code{ff_matrix} object that is used as a block, and is useful for
    distributed computing where each block is processed on a different node.
}
\usage{
\method{as.symDMatrix}{character}(x, ...)
}
\arguments{
    \item{x}{
        A character vector with path names to \code{RData} files.
    }
    \item{...}{
        Additional arguments (currently unused).
    }
}
\details{
    The \code{RData} files must be ordered by block: \code{G11, G12, G13, ...,
    G1q, G22, G23, ..., G2q, ..., Gqq}. The matrix-like objects are
    initialized similarly to \code{load.symDMatrix}.
}
\value{
    A \code{symDMatrix} object.
}
\seealso{
    \code{\link[base]{list.files}} to create a character vector of file paths
    that match a certain pattern.
}
