context("symDMatrix")

library(BEDMatrix)

# Prepare dummy data
X <- BEDMatrix(system.file("extdata", "mice.bed", package = "BEDMatrix"))
G2 <- tcrossprod(scale(X[1:100, ]))
G2 <- G2 / mean(diag(G2))

# Prepare dummy symDMatrix
load.symDMatrix(system.file("extdata", "G.RData", package = "symDMatrix"))

test_that("diag", {

    expect_equal(diag(G), diag(G2))

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

    #expect_equal(G[1], G2[1]) # not yet implemented
    #expect_equal(G[1:2], G2[1:2]) # not yet implemented
    #expect_equal(G[2:1], G2[2:1]) # not yet implemented
    #expect_equal(G[c(3, 1)], G2[c(3, 1)]) # not yet implemented
    #expect_equal(G[c(TRUE, FALSE, TRUE)], G2[c(TRUE, FALSE, TRUE)]) # not yet implemented
    #expect_equal(G[G2 > 1], G2[G2 > 1]) # not yet implemented
    #expect_equal(typeof(G[1]), "double") # not yet implemented

    expect_equal(G["A048005080_A048005080", ], G2["A048005080_A048005080", ])
    expect_equal(G[, "A048005080_A048005080"], G2[, "A048005080_A048005080"])
    expect_equal(G["A048005080_A048005080", "A048005080_A048005080"], G2["A048005080_A048005080", "A048005080_A048005080"])
    expect_equal(G["A048005080_A048005080", , drop = FALSE], G2["A048005080_A048005080", , drop = FALSE])
    expect_equal(G[, "A048005080_A048005080", drop = FALSE], G2[, "A048005080_A048005080", drop = FALSE])
    expect_equal(G["A048005080_A048005080", "A048005080_A048005080", drop = FALSE], G2["A048005080_A048005080", "A048005080_A048005080", drop = FALSE])
    expect_equal(typeof(G["A048005080_A048005080", ]), "double")

    expect_equal(G[c("A048005080_A048005080", "A048006063_A048006063"), ], G2[c("A048005080_A048005080", "A048006063_A048006063"), ])
    expect_equal(G[, c("A048005080_A048005080", "A048006063_A048006063")], G2[, c("A048005080_A048005080", "A048006063_A048006063")])
    expect_equal(G[c("A048005080_A048005080", "A048006063_A048006063"), c("A048005080_A048005080", "A048006063_A048006063")], G2[c("A048005080_A048005080", "A048006063_A048006063"), c("A048005080_A048005080", "A048006063_A048006063")])
    expect_equal(G[c("A048005080_A048005080", "A048006063_A048006063"), , drop = FALSE], G2[c("A048005080_A048005080", "A048006063_A048006063"), , drop = FALSE])
    expect_equal(G[, c("A048005080_A048005080", "A048006063_A048006063"), drop = FALSE], G2[, c("A048005080_A048005080", "A048006063_A048006063"), drop = FALSE])
    expect_equal(G[c("A048005080_A048005080", "A048006063_A048006063"), c("A048005080_A048005080", "A048006063_A048006063"), drop = FALSE], G2[c("A048005080_A048005080", "A048006063_A048006063"), c("A048005080_A048005080", "A048006063_A048006063"), drop = FALSE])
    expect_equal(typeof(G[c("A048005080_A048005080", "A048006063_A048006063"), ]), "double")

    expect_equal(G[c("A048006063_A048006063", "A048005080_A048005080"), ], G2[c("A048006063_A048006063", "A048005080_A048005080"), ])
    expect_equal(G[, c("A048006063_A048006063", "A048005080_A048005080")], G2[, c("A048006063_A048006063", "A048005080_A048005080")])
    expect_equal(G[c("A048006063_A048006063", "A048005080_A048005080"), c("A048006063_A048006063", "A048005080_A048005080")], G2[c("A048006063_A048006063", "A048005080_A048005080"), c("A048006063_A048006063", "A048005080_A048005080")])
    expect_equal(G[c("A048006063_A048006063", "A048005080_A048005080"), , drop = FALSE], G2[c("A048006063_A048006063", "A048005080_A048005080"), , drop = FALSE])
    expect_equal(G[, c("A048006063_A048006063", "A048005080_A048005080"), drop = FALSE], G2[, c("A048006063_A048006063", "A048005080_A048005080"), drop = FALSE])
    expect_equal(G[c("A048006063_A048006063", "A048005080_A048005080"), c("A048006063_A048006063", "A048005080_A048005080"), drop = FALSE], G2[c("A048006063_A048006063", "A048005080_A048005080"), c("A048006063_A048006063", "A048005080_A048005080"), drop = FALSE])
    expect_equal(typeof(G[c("A048006063_A048006063", "A048005080_A048005080"), ]), "double")

    expect_equal(G[c("A048006555_A048006555", "A048005080_A048005080"), ], G2[c("A048006555_A048006555", "A048005080_A048005080"), ])
    expect_equal(G[, c("A048006555_A048006555", "A048005080_A048005080")], G2[, c("A048006555_A048006555", "A048005080_A048005080")])
    expect_equal(G[c("A048006555_A048006555", "A048005080_A048005080"), c("A048006555_A048006555", "A048005080_A048005080")], G2[c("A048006555_A048006555", "A048005080_A048005080"), c("A048006555_A048006555", "A048005080_A048005080")])
    expect_equal(G[c("A048006555_A048006555", "A048005080_A048005080"), , drop = FALSE], G2[c("A048006555_A048006555", "A048005080_A048005080"), , drop = FALSE])
    expect_equal(G[, c("A048006555_A048006555", "A048005080_A048005080"), drop = FALSE], G2[, c("A048006555_A048006555", "A048005080_A048005080"), drop = FALSE])
    expect_equal(G[c("A048006555_A048006555", "A048005080_A048005080"), c("A048006555_A048006555", "A048005080_A048005080"), drop = FALSE], G2[c("A048006555_A048006555", "A048005080_A048005080"), c("A048006555_A048006555", "A048005080_A048005080"), drop = FALSE])
    expect_equal(typeof(G[c("A048006555_A048006555", "A048005080_A048005080"), ]), "double")

})
