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
#' @seealso [symDMatrix-class] for the `symDMatrix` class.
#' @docType package
#' @name symDMatrix-package
#' @aliases symDMatrix-package
#' @import methods
NULL
