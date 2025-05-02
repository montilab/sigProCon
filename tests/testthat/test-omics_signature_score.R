test_that("eset is an expressionset object", {
  eset <- eset
  expect_true(is(eset, "ExpressionSet"))
})

test_that("signature is a list", {
  expect_true(is(signature, "list"))
})

test_that("all signature genes are present in featureNames(eset)") {
  signature <- signature
  expect_true(all(signature[[1]] %in% featureNames(eset)))
}
