source("setup.R")

expect_equal(G2, as.symDMatrix(G2, blockSize = 17, folderOut = testDir())[])
expect_equal(G2, as.symDMatrix(list.files(system.file("extdata", package = "symDMatrix"), pattern = "data_.*\\.RData", full.names = TRUE))[])
