library(BGData)

# Load example genotypes
X <- BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix"))

# Convert BEDMatrix to matrix
X <- as.matrix(X)

# Compute centers and scales
centers <- apply(X, 2, mean, na.rm = TRUE)
scales <- apply(X, 2, sd, na.rm = TRUE)

# Compute G matrix
G1 <- getG(X, center = centers, scale = scales)

# Create symDMatrix object from G matrix
G2 <- as.symDMatrix(G1, folderOut = "fromMatrix")
attr(G2, "centers") <- centers
attr(G2, "scales") <- scales
