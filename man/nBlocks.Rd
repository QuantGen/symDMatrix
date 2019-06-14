\name{nBlocks}
\alias{nBlocks}
\title{Return the Number of Column/Row Blocks of a symDMatrix Object}
\description{
This function returns the number of row blocks the original matrix has been
partitioned into.
}
\usage{
nBlocks(x)
}
\arguments{
\item{x}{A \linkS4class{symDMatrix} object.}
}
\value{
The number of column/row blocks of a \linkS4class{symDMatrix} object.
}
\examples{
# Load example symDMatrix (G)
load.symDMatrix(system.file("extdata", "G.RData", package = "symDMatrix"), readonly = TRUE)

# Get the number of row blocks the original matrix was partitioned into
nBlocks(G)
}
