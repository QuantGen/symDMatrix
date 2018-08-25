# Load example symDMatrix (G)
load.symDMatrix(system.file("extdata", "G.RData", package = "symDMatrix"), readonly = TRUE)

# Get the number of row blocks the original matrix was partitioned into
nBlocks(G)
