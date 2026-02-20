## code to prepare `signature` dataset goes here

set.seed(7)

#create a signature of 100 random genes that are in our eset.
#eset was created in eset.R, and contains 10,000 'genes' names Gene1 ... Gene10000
n_genes <- 10000
gene_names <- paste0("Gene", 1:n_genes)

#take random sample of 100
sample_genes <- sample(gene_names, 100, replace = FALSE)
#make into named list
signature <- list('signature1' = sample_genes)

usethis::use_data(signature, overwrite = TRUE)
