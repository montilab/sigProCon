library(SummarizedExperiment)

test_that("eset is an expressionset object", {
  eset <- ExpressionSet(matrix(rnorm(20), nrow = 5))
  se <- SummarizedExperiment(assays = list(counts = matrix(rnorm(20), nrow = 5)))
  expect_true(inherits(eset, "ExpressionSet") || inherits(eset, "SummarizedExperiment"))
  expect_true(inherits(se, "ExpressionSet") || inherits(se, "SummarizedExperiment"))
})

#this is not the right kind of thing to do
#test_that("signature is a list", {
#  expect_true(is(signature, "list"))
#})

test_that("sig_score length matches ncol(eset)", {
  eset <- ExpressionSet(matrix(rnorm(20), nrow = 5))
  signature <- list('sig1' = sample(rownames(eset), 2, replace = FALSE))
  # Correct length: should not error
  sig_score_good <- setNames(runif(ncol(eset)), sampleNames(eset))
  expect_silent(
    omics_signature_projection(eset, signature, sig_score = sig_score_good)
  )

  # Incorrect length: should error
  sig_score_bad <- setNames(runif(ncol(eset) + 1), c(sampleNames(eset), "extra"))
  expect_error(
    omics_signature_projection(eset, signature, sig_score = sig_score_bad),
    "length\\(sig_score\\)==ncol\\(eset\\)"
  )
})

test_that("omics_signature_projection returns correct output structure", {
  eset <- data(eset)
  signature <- data(signature)
  result <- omics_signature_projection(eset, signature)

  expect_type(result, "list")
  expect_true(all(c("score_cor", "sig_score", "heatmap_all_genes", "heatmap_sig_genes", "ks") %in% names(result)))
  expect_s3_class(result$score_cor, "data.frame")
  expect_s3_class(result$sig_score, "data.frame")
  # Check that score_cor and pval_cor columns exist and are numeric
  expect_true(is.numeric(result$score_cor$score_cor))
  expect_true(is.numeric(result$score_cor$pval_cor))
  expect_true(inherits(result$heatmap_all_genes, "Heatmap"))
  expect_true(inherits(result$heatmap_sig_genes, "Heatmap"))
  expect_true(is.list(result$ks))
})
