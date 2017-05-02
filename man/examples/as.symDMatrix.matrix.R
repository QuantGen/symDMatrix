library(BGData)

# Load example genotypes
X <- BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix"))

# Compute centers and scales
centers <- chunkedApply(X, 2, mean, na.rm = TRUE)
scales <- chunkedApply(X, 2, sd, na.rm = TRUE)

# Compute G matrix
G <- getG(X, center = centers, scale = scales)

# Create symDMatrix object from G matrix
symG <- as.symDMatrix(G, center = centers, scale = scales, folderOut = "fromMatrix")
