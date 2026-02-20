## code to prepare `toy_signature` dataset goes here

usethis::use_data(toy_signature, overwrite = TRUE)

set.seed(123)  # for reproducible results
ngenes <- 100  # number of genes
sig_len <- 10  # signature length

feature_names <- sprintf("gene_%03d", seq(1, ngenes))
toy_signature <- feature_names[c(seq(1, sig_len), sample(seq(sig_len + 1, ngenes), size = 3))]

# Write to to rda and save in the data folder
usethis::use_data(toy_signature, overwrite = TRUE)
