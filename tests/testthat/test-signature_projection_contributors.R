library(SummarizedExperiment)

build_test_eset <- function(n_genes = 30, n_samples = 8) {
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
  eset <- build_test_eset(n_genes = 20, n_samples = 6)
  signature <- list(sig1 = sample(rownames(eset), 5, replace = FALSE))
  # Correct length: should not error
  sig_score_good <- setNames(runif(ncol(eset)), sampleNames(eset))
  expect_silent(
    sigProCon::signature_projection_contributors(eset, signature, sig_score = sig_score_good)
  )

  # Incorrect length: should error
  sig_score_bad <- setNames(runif(ncol(eset) + 1), c(sampleNames(eset), "extra"))
  expect_error(
    sigProCon::signature_projection_contributors(eset, signature, sig_score = sig_score_bad),
    "length\\(sig_score\\) == ncol\\(eset\\)"
  )
})

test_that("signature_projection_contributors returns correct output structure", {
  eset <- build_test_eset(n_genes = 40, n_samples = 10)
  signature <- list(sig1 = sample(rownames(eset), 12, replace = FALSE))
  result <- sigProCon::signature_projection_contributors(eset, signature)

  expect_type(result, "list")
  expect_true(all(c("score_cor", "sig_score", "heatmap_all_genes", "heatmap_sig_genes", "ks") %in% names(result)))
  expect_s3_class(result$score_cor, "data.frame")
  expect_s3_class(result$sig_score, "data.frame")
  # Check that score_cor and pval_cor columns exist and are numeric
  expect_true(is.numeric(result$score_cor$score_cor))
  expect_true(is.numeric(result$score_cor$pval_cor))
  expect_true(any(inherits(result$heatmap_all_genes, c("Heatmap", "HeatmapList"))))
  expect_true(any(inherits(result$heatmap_sig_genes, c("Heatmap", "HeatmapList"))))
  expect_true(is.list(result$ks))
})

test_that("signature_projection_contributors supports eigengene score method", {
  skip_if_wgcna_unavailable()
  eset <- build_test_eset(n_genes = 40, n_samples = 10)
  signature <- list(sig1 = sample(rownames(eset), 12, replace = FALSE))

  result <- sigProCon::signature_projection_contributors(
    eset = eset,
    signature = signature,
    method = "eigengene"
  )

  expect_type(result, "list")
  expect_true(all(c("score_cor", "sig_score", "heatmap_all_genes", "heatmap_sig_genes", "ks") %in% names(result)))
  expect_equal(nrow(result$sig_score), unname(ncol(eset)))
})

test_that("signature_projection_contributors supports pc score method", {
  eset <- build_test_eset(n_genes = 40, n_samples = 10)
  signature <- list(sig1 = sample(rownames(eset), 12, replace = FALSE))

  result <- sigProCon::signature_projection_contributors(
    eset = eset,
    signature = signature,
    method = "pc"
  )

  expect_type(result, "list")
  expect_true(all(c("score_cor", "sig_score", "heatmap_all_genes", "heatmap_sig_genes", "ks") %in% names(result)))
  expect_equal(nrow(result$sig_score), unname(ncol(eset)))
})

test_that("signature_projection_contributors can skip heatmap generation", {
  eset <- build_test_eset(n_genes = 40, n_samples = 10)
  signature <- list(sig1 = sample(rownames(eset), 12, replace = FALSE))

  result <- sigProCon::signature_projection_contributors(
    eset = eset,
    signature = signature,
    make_heatmap_all = FALSE,
    make_heatmap_sig = FALSE
  )

  expect_null(result$heatmap_all_genes)
  expect_null(result$heatmap_sig_genes)
})

test_that("spc_heatmap_* can be called from signature_projection_contributors output", {
  eset <- build_test_eset(n_genes = 40, n_samples = 10)
  signature <- list(sig1 = sample(rownames(eset), 12, replace = FALSE))

  result <- sigProCon::signature_projection_contributors(
    eset = eset,
    signature = signature,
    make_heatmap_all = FALSE,
    make_heatmap_sig = FALSE
  )

  hm_all <- sigProCon::spc_heatmap_all(eset = eset, spc_out = result)
  hm_sig <- sigProCon::spc_heatmap_sig(eset = eset, spc_out = result)

  expect_true(any(inherits(hm_all, c("Heatmap", "HeatmapList"))))
  expect_true(any(inherits(hm_sig, c("Heatmap", "HeatmapList"))))
})

test_that("spc_heatmap_* support additional Heatmap arguments", {
  eset <- build_test_eset(n_genes = 40, n_samples = 10)
  signature <- list(sig1 = sample(rownames(eset), 12, replace = FALSE))

  result <- sigProCon::signature_projection_contributors(
    eset = eset,
    signature = signature,
    make_heatmap_all = FALSE,
    make_heatmap_sig = FALSE
  )

  hm_all <- sigProCon::spc_heatmap_all(
    eset = eset,
    spc_out = result,
    row_names_side = "right",
    show_column_names = TRUE
  )
  hm_sig <- sigProCon::spc_heatmap_sig(
    eset = eset,
    spc_out = result,
    row_names_side = "right",
    show_column_names = TRUE
  )

  expect_true(any(inherits(hm_all, c("Heatmap", "HeatmapList"))))
  expect_true(any(inherits(hm_sig, c("Heatmap", "HeatmapList"))))
})

test_that("spc_heatmap_all subsample keeps top MAD genes plus signature genes", {
  set.seed(123)
  mat <- matrix(
    c(
      rep(1, 8), # gene1: low MAD, in signature
      1:8,       # gene2: high MAD
      rep(c(0, 10), 4), # gene3: high MAD
      c(1, 2, 1, 2, 1, 2, 1, 2), # gene4: moderate MAD
      c(5, 5, 5, 5, 5, 5, 6, 6), # gene5: low MAD
      c(2, 2, 3, 3, 4, 4, 5, 5)  # gene6: moderate MAD
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(paste0("gene", 1:6), paste0("sample", 1:8))
  )
  eset <- Biobase::ExpressionSet(mat)
  signature <- list(sig1 = "gene1")

  result <- sigProCon::signature_projection_contributors(
    eset = eset,
    signature = signature,
    make_heatmap_all = FALSE,
    make_heatmap_sig = FALSE
  )

  prep <- sigProCon:::.spc_prepare_heatmap_data(
    eset = eset,
    spc_out = result,
    subsample = 2
  )

  selected_genes <- rownames(Biobase::exprs(prep$eset_srt))
  expect_true("gene1" %in% selected_genes)
  expect_true(length(selected_genes) >= 2)
  expect_true(length(selected_genes) <= 3)
})
