build_score_test_eset <- function(n_genes = 30, n_samples = 8) {
  set.seed(1)
  mat <- matrix(
    rnorm(n_genes * n_samples),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(paste0("gene", seq_len(n_genes)), paste0("sample", seq_len(n_samples)))
  )
  Biobase::ExpressionSet(mat)
}

skip_if_wgcna_unavailable <- function() {
  rscript <- Sys.which("Rscript")
  if (nzchar(rscript)) {
    tf <- tempfile(fileext = ".R")
    on.exit(unlink(tf), add = TRUE)
    writeLines("suppressPackageStartupMessages(library(WGCNA));cat('ok')", tf)
    status <- system2(
      rscript,
      tf,
      stdout = FALSE,
      stderr = FALSE
    )
    if (!identical(status, 0L)) {
      testthat::skip("WGCNA runtime is unavailable in this environment")
    }
  } else {
    testthat::skip_if_not_installed("WGCNA")
  }
}

test_that("omics_signature_score returns expected score vector for GSVA", {
  eset <- build_score_test_eset(n_genes = 30, n_samples = 6)
  signature <- list(sig1 = paste0("gene", 1:8))

  score_gsva <- sigProCon::omics_signature_score(
    eset = eset,
    signature = signature,
    method = "GSVA"
  )

  expect_type(score_gsva, "double")
  expect_length(score_gsva, ncol(eset))
  expect_identical(names(score_gsva), Biobase::sampleNames(eset))
})

test_that("omics_signature_score returns expected score vector for pc", {
  eset <- build_score_test_eset(n_genes = 30, n_samples = 6)
  signature <- list(sig1 = paste0("gene", 1:8))

  score_pc <- sigProCon::omics_signature_score(
    eset = eset,
    signature = signature,
    method = "pc"
  )

  expect_type(score_pc, "double")
  expect_length(score_pc, ncol(eset))
  expect_identical(names(score_pc), Biobase::sampleNames(eset))
})

test_that("omics_signature_score returns expected score vector for eigengene", {
  skip_if_wgcna_unavailable()
  eset <- build_score_test_eset(n_genes = 30, n_samples = 6)
  signature <- list(sig1 = paste0("gene", 1:8))

  score_eigengene <- sigProCon::omics_signature_score(
    eset = eset,
    signature = signature,
    method = "eigengene"
  )

  expect_type(score_eigengene, "double")
  expect_length(score_eigengene, ncol(eset))
  expect_identical(names(score_eigengene), Biobase::sampleNames(eset))
})

test_that("omics_signature_score validates signature features", {
  eset <- build_score_test_eset(n_genes = 20, n_samples = 5)
  signature <- list(sig1 = c("gene1", "gene2", "missing_gene"))

  expect_error(
    sigProCon::omics_signature_score(eset = eset, signature = signature, method = "eigengene")
  )
})

test_that("extract_eigengenes works with a single signature and method multiple", {
  skip_if_wgcna_unavailable()
  eset <- build_score_test_eset(n_genes = 25, n_samples = 7)
  signature <- list(sig1 = paste0("gene", 1:6))
  eigenset <- sigProCon::extract_eigengenes(eset = eset, signatures = signature, method = "multiple")

  expect_true(methods::is(eigenset, "ExpressionSet"))
  expect_identical(Biobase::sampleNames(eigenset), Biobase::sampleNames(eset))
  expect_identical(Biobase::featureNames(eigenset), "sig1")
})
