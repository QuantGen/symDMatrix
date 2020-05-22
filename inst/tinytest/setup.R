testDir <- function() {
    paste0(tempdir(), "/symDMatrix-", symDMatrix:::randomString(), "/")
}

# Prepare dummy data
X <- suppressMessages(BEDMatrix::BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix")))
X <- scale(X)
X[is.na(X)] <- 0
G2 <- tcrossprod(X)
G2 <- G2 / mean(diag(G2))

# Prepare dummy symDMatrix
suppressMessages(load.symDMatrix(system.file("extdata", "G.RData", package = "symDMatrix"), readonly = TRUE))
