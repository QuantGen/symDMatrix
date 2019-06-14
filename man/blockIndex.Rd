\name{blockIndex}
\alias{blockIndex}
\title{Return the Block Structure of a symDMatrix Object}
\description{
This function returns the block structure of a \code{symDMatrix} object and
can be useful when implementing custom indexing techniques.
}
\usage{
blockIndex(x)
}
\arguments{
  \item{x}{
    A \code{symDMatrix} object.
  }
}
\value{
A matrix with three columns: the block number, the start index and the end
index.
}
