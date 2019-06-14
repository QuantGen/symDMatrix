\docType{package}
\name{symDMatrix-package}
\alias{symDMatrix-package}
\title{A Package Providing Symmetric Matrices Partitioned into File-Backed Blocks.}
\description{
A Package Providing Symmetric Matrices Partitioned into File-Backed Blocks.
}
\section{Example Dataset}{
 The example dataset in the \code{extdata} folder is the
G matrix of the dummy dataset that comes with the
\link[BEDMatrix:BEDMatrix-package]{BEDMatrix} package. It has been generated as
follows:\preformatted{library(BGData)
X <- BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix"))
G <- getG_symDMatrix(X, blockSize = 17, folderOut = "inst/extdata")
}

To load the dataset:\preformatted{load.symDMatrix(system.file("extdata", "G.RData", package = "symDMatrix"),
readonly = TRUE)
}

To demonstrate the \code{\link[=as.symDMatrix]{as.symDMatrix()}} method for character vectors, \code{RData}
files for each block have been generated:\preformatted{for (i in 1:nBlocks(G)) {
    for (j in i:nBlocks(G)) {
        block <- G[[i]][[j]]
        save(block, file = paste0("inst/extdata/data_", i, "_", j, ".RData"))
    }
}
}
}

\seealso{
\linkS4class{symDMatrix} for the \code{symDMatrix} class.
}
