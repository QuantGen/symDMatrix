source("setup.R")

expect_equal(blockSize(G), 17)
expect_equal(blockSize(G, last = FALSE), 17)
expect_equal(blockSize(G, last = TRUE), 16)
