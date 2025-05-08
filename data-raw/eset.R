## code to prepare `eset` dataset goes here

usethis::use_data(eset, overwrite = TRUE)

set.seed(7)

# Simulate expression matrix (genes x samples)
n_genes <- 10000
n_samples <- 30

exprs <- matrix(
  rpois(n_genes * n_samples, lambda = 100),  # random counts
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(
    paste0("Gene", 1:n_genes),
    paste0("Sample", 1:n_samples)
  )
)

# Make normalized count matrix
# Calculate library sizes (total counts per sample)
lib_sizes <- colSums(exprs)
cpm <- t(t(exprs) / lib_sizes * 1e6)
logcpm <- log2(cpm + 1)

# Simulate sample metadata
phenoData <- data.frame(
  Group = rep(c("Control", "Treated"), each = 15),
  row.names = paste0("Sample", 1:n_samples)
)

# Create ExpressionSet
eset <- ExpressionSet(
  assayData = logcpm,
  phenoData = AnnotatedDataFrame(phenoData)
)


# Write to to rda and save in the data folder
usethis::use_data(eset, overwrite = TRUE)
