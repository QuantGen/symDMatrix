# Load example symDMatrix (G)
load.symDMatrix(system.file("extdata", "G.RData", package = "symDMatrix"))

# Get the block size
blockSize(G)

# Get the block size of the trailing blocks
blockSize(G, last = TRUE)
