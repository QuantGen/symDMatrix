\name{load.symDMatrix}
\alias{load.symDMatrix}
\title{Load symDMatrix Objects from .RData Files}
\description{
This function is similar to \code{load}, but it also initializes the
\code{ff_matrix} blocks in the \code{symDMatrix} object.
}
\usage{
load.symDMatrix(file, readonly = FALSE, envir = parent.frame())
}
\arguments{
  \item{file}{
    The name of an .RData file to be loaded.
  }
  \item{readonly}{
    Set to TRUE to forbid writing to existing files.
  }
  \item{envir}{
    The environment where to load the data.
  }
}
