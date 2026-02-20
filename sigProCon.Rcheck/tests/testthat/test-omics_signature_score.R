test_that("eset is an expressionset object", {
  eset <- Biobase::ExpressionSet(
    matrix(rnorm(40), nrow = 10, dimnames = list(paste0("gene", 1:10), paste0("sample", 1:4)))
  )
  expect_true(is(eset, "ExpressionSet"))
})

test_that("signature is a list", {
  signature <- list(signature1 = c("gene1", "gene2", "gene3"))
  expect_true(is(signature, "list"))
})

test_that("all signature genes are present in featureNames(eset)", {
  eset <- Biobase::ExpressionSet(
    matrix(rnorm(40), nrow = 10, dimnames = list(paste0("gene", 1:10), paste0("sample", 1:4)))
  )
  signature <- list(signature1 = c("gene1", "gene2", "gene3"))
  expect_true(all(signature[[1]] %in% featureNames(eset)))
})
