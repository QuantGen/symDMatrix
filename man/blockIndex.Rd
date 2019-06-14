\name{blockIndex}
\alias{blockIndex}
\title{Return the Block Structure of a symDMatrix Object}
\usage{
blockIndex(x)
}
\arguments{
\item{x}{A \linkS4class{symDMatrix} object.}
}
\value{
A matrix with three columns: the block number, the start index and
the end index.
}
\description{
This function returns the block structure of a \linkS4class{symDMatrix} object and
can be useful when implementing custom indexing techniques.
}
