context("symDMatrix")

# Prepare dummy data
X <- suppressMessages(BEDMatrix::BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix")))
X <- scale(X)
X[is.na(X)] <- 0
G2 <- tcrossprod(X)
G2 <- G2 / mean(diag(G2))

# Prepare dummy symDMatrix
suppressMessages(load.symDMatrix(system.file("extdata", "G.RData", package = "symDMatrix")))

test_that("symDMatrix", {

    # Test that there is at least one block
    expect_error(symDMatrix(data = list()), "data needs to contain at least one block")

    # Test that data has the right structure
    expect_error(symDMatrix(data = list(list(), list(), list())), "data needs to be a nested list in the following structure")

    # Test that all blocks are matrix-like objects
    listBlock <- vector(mode = "list", length = 5)
    expect_error(symDMatrix(data = list(list(listBlock, listBlock, listBlock), list(listBlock, listBlock), list(listBlock))), "data: all blocks need to be matrix-like objects")

    # Test that all blocks per row have the same number of rows
    ffBlock <- ff::ff(initdata = rnorm(25), dim = c(5, 5))
    wrongRows <- ff::ff(initdata = rnorm(15), dim = c(3, 5))
    expect_error(symDMatrix(data = list(list(ffBlock, wrongRows, ffBlock), list(ffBlock, ffBlock), list(ffBlock))), "data: all blocks per row need the same number of rows")

    # Test that all blocks per column have the same number of columns
    wrongColumns <- ff::ff(initdata = rnorm(15), dim = c(5, 3))
    expect_error(symDMatrix(data = list(list(ffBlock, wrongColumns, ffBlock), list(ffBlock, ffBlock), list(ffBlock))), "data: all blocks per column need the same number of columns")

    # Test that the first block is square
    notSquare <- ff::ff(initdata = rnorm(20), dim = c(4, 5))
    expect_error(symDMatrix(data = list(list(notSquare))), "data: the first block needs to be square")

    # Test that all non-final blocks need to be square
    expect_error(symDMatrix(data = list(list(ffBlock, wrongColumns, ffBlock), list(wrongColumns, ffBlock), list(ffBlock))), "data: non-final blocks need to be square")

})

test_that("diag", {

    expect_equal(diag(G), diag(G2))

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

test_that("subsetting", {

    expect_equal(G[], G2[])
    expect_equal(typeof(G[]), "double")

    expect_equal(G[1, ], G2[1, ])
    expect_equal(G[, 1], G2[, 1])
    expect_equal(G[1, 1], G2[1, 1])
    expect_equal(G[1, , drop = FALSE], G2[1, , drop = FALSE])
    expect_equal(G[, 1, drop = FALSE], G2[, 1, drop = FALSE])
    expect_equal(G[1, 1, drop = FALSE], G2[1, 1, drop = FALSE])
    expect_equal(typeof(G[1, ]), "double")

    expect_equal(G[1:2, ], G2[1:2, ])
    expect_equal(G[, 1:2], G2[, 1:2])
    expect_equal(G[1:2, 1:2], G2[1:2, 1:2])
    expect_equal(G[1:2, , drop = FALSE], G2[1:2, , drop = FALSE])
    expect_equal(G[, 1:2, drop = FALSE], G2[, 1:2, drop = FALSE])
    expect_equal(G[1:2, 1:2, drop = FALSE], G2[1:2, 1:2, drop = FALSE])
    expect_equal(typeof(G[1:2, ]), "double")

    expect_equal(G[2:1, ], G2[2:1, ])
    expect_equal(G[, 2:1], G2[, 2:1])
    expect_equal(G[2:1, 2:1], G2[2:1, 2:1])
    expect_equal(G[2:1, , drop = FALSE], G2[2:1, , drop = FALSE])
    expect_equal(G[, 2:1, drop = FALSE], G2[, 2:1, drop = FALSE])
    expect_equal(G[2:1, 2:1, drop = FALSE], G2[2:1, 2:1, drop = FALSE])
    expect_equal(typeof(G[2:1, ]), "double")

    expect_equal(G[c(3, 1), ], G2[c(3, 1), ])
    expect_equal(G[, c(3, 1)], G2[, c(3, 1)])
    expect_equal(G[c(3, 1), c(3, 1)], G2[c(3, 1), c(3, 1)])
    expect_equal(G[c(3, 1), , drop = FALSE], G2[c(3, 1), , drop = FALSE])
    expect_equal(G[, c(3, 1), drop = FALSE], G2[, c(3, 1), drop = FALSE])
    expect_equal(G[c(3, 1), c(3, 1), drop = FALSE], G2[c(3, 1), c(3, 1), drop = FALSE])
    expect_equal(typeof(G[c(3, 1), ]), "double")

    expect_equal(G[c(TRUE, FALSE, TRUE), ], G2[c(TRUE, FALSE, TRUE), ])
    expect_equal(G[, c(TRUE, FALSE, TRUE)], G2[, c(TRUE, FALSE, TRUE)])
    expect_equal(G[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)], G2[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)])
    expect_equal(G[c(TRUE, FALSE, TRUE), , drop = FALSE], G2[c(TRUE, FALSE, TRUE), , drop = FALSE])
    expect_equal(G[, c(TRUE, FALSE, TRUE), drop = FALSE], G2[, c(TRUE, FALSE, TRUE), drop = FALSE])
    expect_equal(G[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE], G2[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE])
    expect_equal(typeof(G[c(TRUE, FALSE, TRUE), ]), "double")

    expect_equal(G[1], G2[1])
    expect_equal(G[1:2], G2[1:2])
    expect_equal(G[2:1], G2[2:1])
    expect_equal(G[c(3, 1)], G2[c(3, 1)])
    #expect_equal(G[c(TRUE, FALSE, TRUE)], G2[c(TRUE, FALSE, TRUE)]) # not yet implemented
    #expect_equal(G[G2 > 1], G2[G2 > 1]) # not yet implemented
    expect_equal(typeof(G[1]), "double")

    expect_equal(G["per0_per0", ], G2["per0_per0", ])
    expect_equal(G[, "per0_per0"], G2[, "per0_per0"])
    expect_equal(G["per0_per0", "per0_per0"], G2["per0_per0", "per0_per0"])
    expect_equal(G["per0_per0", , drop = FALSE], G2["per0_per0", , drop = FALSE])
    expect_equal(G[, "per0_per0", drop = FALSE], G2[, "per0_per0", drop = FALSE])
    expect_equal(G["per0_per0", "per0_per0", drop = FALSE], G2["per0_per0", "per0_per0", drop = FALSE])
    expect_equal(typeof(G["per0_per0", ]), "double")

    expect_equal(G[c("per0_per0", "per1_per1"), ], G2[c("per0_per0", "per1_per1"), ])
    expect_equal(G[, c("per0_per0", "per1_per1")], G2[, c("per0_per0", "per1_per1")])
    expect_equal(G[c("per0_per0", "per1_per1"), c("per0_per0", "per1_per1")], G2[c("per0_per0", "per1_per1"), c("per0_per0", "per1_per1")])
    expect_equal(G[c("per0_per0", "per1_per1"), , drop = FALSE], G2[c("per0_per0", "per1_per1"), , drop = FALSE])
    expect_equal(G[, c("per0_per0", "per1_per1"), drop = FALSE], G2[, c("per0_per0", "per1_per1"), drop = FALSE])
    expect_equal(G[c("per0_per0", "per1_per1"), c("per0_per0", "per1_per1"), drop = FALSE], G2[c("per0_per0", "per1_per1"), c("per0_per0", "per1_per1"), drop = FALSE])
    expect_equal(typeof(G[c("per0_per0", "per1_per1"), ]), "double")

    expect_equal(G[c("per1_per1", "per0_per0"), ], G2[c("per1_per1", "per0_per0"), ])
    expect_equal(G[, c("per1_per1", "per0_per0")], G2[, c("per1_per1", "per0_per0")])
    expect_equal(G[c("per1_per1", "per0_per0"), c("per1_per1", "per0_per0")], G2[c("per1_per1", "per0_per0"), c("per1_per1", "per0_per0")])
    expect_equal(G[c("per1_per1", "per0_per0"), , drop = FALSE], G2[c("per1_per1", "per0_per0"), , drop = FALSE])
    expect_equal(G[, c("per1_per1", "per0_per0"), drop = FALSE], G2[, c("per1_per1", "per0_per0"), drop = FALSE])
    expect_equal(G[c("per1_per1", "per0_per0"), c("per1_per1", "per0_per0"), drop = FALSE], G2[c("per1_per1", "per0_per0"), c("per1_per1", "per0_per0"), drop = FALSE])
    expect_equal(typeof(G[c("per1_per1", "per0_per0"), ]), "double")

    expect_equal(G[c("per3_per3", "per0_per0"), ], G2[c("per3_per3", "per0_per0"), ])
    expect_equal(G[, c("per3_per3", "per0_per0")], G2[, c("per3_per3", "per0_per0")])
    expect_equal(G[c("per3_per3", "per0_per0"), c("per3_per3", "per0_per0")], G2[c("per3_per3", "per0_per0"), c("per3_per3", "per0_per0")])
    expect_equal(G[c("per3_per3", "per0_per0"), , drop = FALSE], G2[c("per3_per3", "per0_per0"), , drop = FALSE])
    expect_equal(G[, c("per3_per3", "per0_per0"), drop = FALSE], G2[, c("per3_per3", "per0_per0"), drop = FALSE])
    expect_equal(G[c("per3_per3", "per0_per0"), c("per3_per3", "per0_per0"), drop = FALSE], G2[c("per3_per3", "per0_per0"), c("per3_per3", "per0_per0"), drop = FALSE])
    expect_equal(typeof(G[c("per3_per3", "per0_per0"), ]), "double")

})
