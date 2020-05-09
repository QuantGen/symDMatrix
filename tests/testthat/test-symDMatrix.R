context("symDMatrix")

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

test_that("symDMatrix", {

    # Test that there is at least one block (will be detected by LinkedMatrix)
    expect_error(symDMatrix())

    # Test that blocks are of type ColumnLinkedMatrix
    matrixSquareBlock <- matrix(data = rnorm(25), nrow = 5, ncol = 5)
    expect_error(symDMatrix(matrixSquareBlock), "blocks need to be of type ColumnLinkedMatrix")

    # Test that number of blocks is consistent (will be detected by LinkedMatrix)
    expect_error(symDMatrix(LinkedMatrix::ColumnLinkedMatrix(matrixSquareBlock, matrixSquareBlock), LinkedMatrix::ColumnLinkedMatrix(matrixSquareBlock)))

    # Test that nested blocks inherit from ff_matrix
    expect_error(symDMatrix(LinkedMatrix::ColumnLinkedMatrix(matrixSquareBlock)), "nested blocks need to inherit from ff_matrix")

    # Test that all blocks per row have the same number of rows
    ffSquareBlock <- ff::ff(dim = c(5, 5), initdata = rnorm(25))
    ffWrongRows <- ff::ff(dim = c(3, 5), initdata = rnorm(15))
    expect_error(symDMatrix(LinkedMatrix::ColumnLinkedMatrix(ffWrongRows, ffWrongRows), LinkedMatrix::ColumnLinkedMatrix(ffSquareBlock, ffSquareBlock)), "non-final blocks need to be square")

    # Test that all blocks per column have the same number of columns (will be detected by LinkedMatrix)
    ffWrongColumns <- ff::ff(dim = c(5, 3), initdata = rnorm(15))
    expect_error(symDMatrix(LinkedMatrix::ColumnLinkedMatrix(ffWrongColumns, ffWrongColumns), LinkedMatrix::ColumnLinkedMatrix(ffSquareBlock, ffSquareBlock)))

    # Test that the first block is square
    ffNotSquare <- ff::ff(dim = c(4, 5), initdata = rnorm(20))
    expect_error(symDMatrix(LinkedMatrix::ColumnLinkedMatrix(ffNotSquare)), "the first block needs to be square")

    # Test that all non-final blocks need to be square
    expect_error(symDMatrix(LinkedMatrix::ColumnLinkedMatrix(ffSquareBlock, ffWrongColumns, ffSquareBlock), LinkedMatrix::ColumnLinkedMatrix(ffSquareBlock, ffWrongColumns, ffSquareBlock), LinkedMatrix::ColumnLinkedMatrix(ffSquareBlock, ffWrongColumns, ffSquareBlock)), "non-final blocks need to be square")

    # Test that matrices are the same
    expect_equal(G2, G[])

})

test_that("diag", {

    expect_equal(diag(G2), diag(G))

})

test_that("nBlocks", {

    expect_equal(nBlocks(G), 3)

})

test_that("blockSize", {

    expect_equal(blockSize(G), 17)
    expect_equal(blockSize(G, last = FALSE), 17)
    expect_equal(blockSize(G, last = TRUE), 16)

})

test_that("blockIndex", {

    expect_equal(blockIndex(G), matrix(data = c(1, 1, 17, 2, 18, 34, 3, 35, 50), ncol = 3, byrow = TRUE, dimnames = list(NULL, c("block", "ini", "end"))))

})

test_that("as.symDMatrix", {

    expect_equal(G2, as.symDMatrix(G2, blockSize = 17, folderOut = testDir())[])
    expect_equal(G2, as.symDMatrix(list.files(system.file("extdata", package = "symDMatrix"), pattern = "data_.*\\.RData", full.names = TRUE))[])

})
