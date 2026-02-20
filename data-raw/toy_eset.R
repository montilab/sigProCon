## code to prepare `toy_eset` dataset goes here

usethis::use_data(toy_eset, overwrite = TRUE)

set.seed(123)  # for reproducible results
nsamples <- 20 # number of samples
ngenes <- 100  # number of genes
sig_len <- 10  # signature length

sig_score <- rnorm(nsamples)

# Simulate sample metadata
phenoData <- data.frame(
  Group = sample(c("head", "tail"), nsamples, replace = TRUE),
  row.names = sprintf("sample_%02d", seq(1, nsamples))
)
## Make expression matrix
dat_matrix <- rbind(
  ## sig_len genes have expression levels correlated to sig_score
  t(replicate(sig_len, rnorm(nsamples, mean = sig_score))),
  ## remaining genes have rnorm(0,1) expression
  matrix(rnorm(n = nsamples * (ngenes - sig_len)), ncol = nsamples),
)
## assign gene and sample names
rownames(dat_matrix) <- sprintf("gene_%03d", seq(1, ngenes))
names(sig_score) <- colnames(dat_matrix) <- sprintf("sample_%02d", seq(1, nsamples))

## create ExpressionSet
toy_eset <- Biobase::ExpressionSet(
  assayData = dat_matrix,
  phenoData = AnnotatedDataFrame(data.frame(sig_score = sig_score, pheno = rnd_pheno))
)
toy_signature <- featureNames(toy_eset)[c(seq(1, sig_len), sample(seq(sig_len + 1, ngenes), size = 3))]

# Write to to rda and save in the data folder
usethis::use_data(toy_eset, overwrite = TRUE)

