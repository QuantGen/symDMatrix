source("setup.R")

expect_equal(blockIndex(G), matrix(data = c(1, 1, 17, 2, 18, 34, 3, 35, 50), ncol = 3, byrow = TRUE, dimnames = list(NULL, c("block", "ini", "end"))))
