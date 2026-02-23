build_pc_test_eset <- function(n_genes = 20, n_samples = 6) {
  set.seed(123)
  mat <- matrix(
    rnorm(n_genes * n_samples),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(paste0("gene", seq_len(n_genes)), paste0("sample", seq_len(n_samples)))
  )
  Biobase::ExpressionSet(mat)
}

test_that("extract_pc returns named vector with sample-level scores", {
  eset <- build_pc_test_eset()
  sig <- paste0("gene", 1:8)

  pc1 <- sigProCon::extract_pc(eset = eset, sig = sig, i = 1)
  manual <- stats::prcomp(t(Biobase::exprs(eset)[sig, , drop = FALSE]), center = TRUE, scale. = FALSE)$x[, 1]

  expect_type(pc1, "double")
  expect_length(pc1, ncol(eset))
  expect_identical(names(pc1), Biobase::sampleNames(eset))
  expect_gt(abs(stats::cor(pc1, manual)), 0.999999)
})

test_that("extract_pc validates inputs", {
  eset <- build_pc_test_eset()

  expect_error(sigProCon::extract_pc(eset = matrix(1:4, 2, 2), sig = "gene1", i = 1))
  expect_error(sigProCon::extract_pc(eset = eset, sig = c("gene1", "missing"), i = 1))
  expect_error(sigProCon::extract_pc(eset = eset, sig = "gene1", i = 0))
  expect_error(sigProCon::extract_pc(eset = eset, sig = "gene1", i = 999))
})
