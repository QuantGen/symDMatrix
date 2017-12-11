#' A Package Providing Symmetric Matrices Assembled from Memory-Mapped Blocks.
#'
#' @section Example Dataset: The example dataset in the `extdata` folder is the
#' G matrix of the dummy dataset that comes with the
#' [BEDMatrix][BEDMatrix::BEDMatrix-package] package. It has been generated as
#' follows:
#'
#' ```
#' library(BGData)
#' X <- BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix"))
#' G <- getG_symDMatrix(X, blockSize = 17, folderOut = "inst/extdata")
#' ```
#'
#' To load the dataset:
#'
#' ```
#' load.symDMatrix(system.file("extdata", "G.RData", package = "symDMatrix"))
#' ```
#'
#' To demonstrate the `as.symDMatrix.character()` method, `RData` files for
#' each block have been generated:
#'
#' ```
#' for (i in 1:nBlocks(G)) {
#'     for (j in i:nBlocks(G)) {
#'         block <- G[[i]][[j]]
#'         save(block, file = paste0("inst/extdata/data_", i, "_", j, ".RData"))
#'     }
#' }
#' ```
#'
#' @seealso [symDMatrix-class] for the `symDMatrix` class.
#' @docType package
#' @name symDMatrix-package
#' @aliases symDMatrix-package
#' @import methods
#' @importClassesFrom LinkedMatrix RowLinkedMatrix
NULL
