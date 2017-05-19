library(BGData)

# Load example genotypes
X <- BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix"))

# Compute centers and scales
centers <- chunkedApply(X, 2, mean, na.rm = TRUE)
scales <- chunkedApply(X, 2, sd, na.rm = TRUE)

# Partition genotypes into three blocks
nBlocks <- 3
blockSize <- ceiling(nrow(X) / nBlocks)

# Create index
i <- 1:nrow(X)
blockIndices <- split(i, ceiling(i / blockSize))

# Iterate block by block
for (r in 1:nBlocks) {
    for (s in r:nBlocks) {
        blockName <- paste0("G_", r, "_", s - r + 1)
        block <- getG(X, center = centers, scale = scales, scaleG = TRUE,
                      i = blockIndices[[r]], i2 = blockIndices[[s]])
        block <- ff::as.ff(block, filename = paste0(blockName, ".bin"), vmode = "double")
        save(block, file = paste0(blockName, ".RData"))
    }
}

# Create symDMatrix from files
G <- as.symDMatrix(list.files(pattern = ".RData$"), centers = centers, scales = scales)
