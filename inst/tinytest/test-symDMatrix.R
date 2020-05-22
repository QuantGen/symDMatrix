source("setup.R")

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
