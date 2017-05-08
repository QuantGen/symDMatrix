# The BGData package is not on CRAN yet and has to be installed from GitHub:
# https://github.com/QuantGen/BGData
if (require("BGData")) {

library(BGData)

# Load example genotypes
X <- BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix"))

# Compute centers and scales
centers <- chunkedApply(X, 2, mean, na.rm = TRUE)
scales <- chunkedApply(X, 2, sd, na.rm = TRUE)

# Compute G matrix
G1 <- getG(X, center = centers, scale = scales)

# Create symDMatrix object from G matrix
G2 <- as.symDMatrix(G1, center = centers, scale = scales, folderOut = "fromMatrix")

}
