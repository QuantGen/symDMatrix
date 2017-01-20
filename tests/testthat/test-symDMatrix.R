context("symDMatrix")

library(BEDMatrix)

# Prepare dummy data
X <- suppressMessages(BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix")))
X <- scale(X)
X[is.na(X)] <- 0
G2 <- tcrossprod(X)
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
