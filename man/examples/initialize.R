# Load example symDMatrix (G)
load.symDMatrix(system.file("extdata", "G.RData", package = "symDMatrix"))

# Create a symDMatrix from a single block
G1 <- symDMatrix(data = list(list(G[, ])))
nBlocks(G1)
blockSize(G1)

# Create a symDMatrix from three blocks (by pretending to partition the
# original matrix into two row blocks)
G2 <- symDMatrix(data = list(list(G[1:25, 1:25], G[1:25, 26:50]), list(G[26:50, 26:50])))
nBlocks(G2)
blockSize(G2)
