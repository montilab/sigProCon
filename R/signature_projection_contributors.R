#' This function generates an annotated heatmap that shows the contributions of genes to the given signature score
#'
#' @param eset an expression set
#' @param signature a list of signatures (at least 1)
#' @param sig_score an aggregate signature score
#' @param cor_method correlation method used by \code{psych::corr.test}
#' @param col_ha a ComplexHeatmap::heatmapAnnotation object with columns' (i.e., samples') annotation
#' @param min_sigsize minimum number of signature genes required to plot signature-only heatmap
#' @param method how to compute the aggregate score (`"GSVA"`, `"eigengene"`, or `"pc"`)
#' @param name name for the heatmap
#' @param gsea logical, whether to use weighted KS statistic
#' @param make_heatmap_all logical, whether to generate all-genes heatmap
#' @param make_heatmap_sig logical, whether to generate signature-only heatmap
#' @param ... additional parameters to pass to ComplexHeatmap::Heatmap
#'
#' @return a list of a dataframe and two heatmaps
#'
#' @import Biobase
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation anno_barplot
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom stats cor
#' @importFrom methods is
#' @importFrom grid gpar
#' @importFrom psych corr.test
#'
#' @export
signature_projection_contributors <- function(
    eset,
    signature,
    sig_score = NULL,
    cor_method = c("pearson", "spearman", "kendall"),
    col_ha = NULL,
    min_sigsize = 3,
    method = c("GSVA", "eigengene", "pc"),
    name = "expression",
    gsea = FALSE,
    make_heatmap_all = FALSE,
    make_heatmap_sig = FALSE,
    ...
) {
  ## BEGIN input checks
  cor_method = match.arg(cor_method)
  stopifnot( methods::is(signature, "list") )
  stopifnot( methods::is(eset, "SummarizedExperiment") || methods::is(eset, "ExpressionSet") )
  stopifnot( is.null(sig_score) || length(sig_score)==ncol(eset) )
  stopifnot( is.null(sig_score) || isTRUE(all.equal(names(sig_score), Biobase::sampleNames(eset))) )
  stopifnot( is.null(col_ha) || isTRUE(all(rownames(col_ha) %in% Biobase::sampleNames(eset))))
  ## END input checks

  if ( methods::is(eset, "SummarizedExperiment") ) {
    eset <- Biobase::ExpressionSet(
      assayData = SummarizedExperiment::assay(eset),
      phenoData = Biobase::AnnotatedDataFrame(colData(eset)),
      featureData = Biobase::AnnotatedDataFrame(rowData(eset))
    )
  }
  ## compute score if not provided
  if ( is.null(sig_score) ) {
    sig_score <- omics_signature_score( eset = eset, signature = signature, method = method)
  }
  ## add signature score to eset metadata
  eset$sig_score <- sig_score
  stopifnot( all(!is.na(eset$sig_score)) )

  ## add correlation (and p-value) of each gene with signature score to fData
  COR <- psych::corr.test(eset$sig_score, t(Biobase::exprs(eset)), method = cor_method)
  stopifnot(nrow(COR$r) == 1)
  stopifnot(nrow(COR$p) == 1)
  Biobase::fData(eset)$score_cor <- drop(COR$r)
  Biobase::fData(eset)$pval_cor <- drop(COR$p)
  Biobase::fData(eset)$insig <- factor(
    ifelse(Biobase::featureNames(eset) %in% signature[[1]], 'signature', 'background'),
    levels = c("signature","background")
  )
  ## prepare sorted object and enrichment output for downstream plotting
  eset_srt <- eset[
    order(Biobase::fData(eset)$score_cor, decreasing = TRUE), # high to low correlation
    order(eset$sig_score, decreasing = TRUE)                  # high to low sig score
  ]
  ks_out <-   .kstest(
    n.x = nrow(eset),
    y = rank(-Biobase::fData(eset)$score_cor)[Biobase::fData(eset)$insig == "signature"],
    weights = if (gsea) Biobase::fData(eset_srt)$score_cor,
    plotting = TRUE
  )
  ## from idx to names
  fData(eset_srt)$leading_edge <- NA
  if (is.null(ks_out$leading_edge)) {
    ks_out$hits <- NA
  } else if (!is.null(ks_out$leading_edge) && ks_out$leading_edge == 0) {
    ks_out$hits <- NA
  } else {
    leading_genes <- Biobase::fData(eset_srt)
    leading_genes <- leading_genes[seq_len(ks_out$leading_edge), , drop = FALSE]
    leading_genes <- leading_genes[leading_genes$insig == "signature", , drop = FALSE]
    ks_out$hits <- rownames(leading_genes)
    ## add information about hits in leading edge
    Biobase::fData(eset_srt)$leading_edge <-
      factor(ifelse(Biobase::featureNames(eset_srt) %in% ks_out$hits, "yes", "no"),
             levels = c("yes", "no"))
  }
  spc_out <- list(
    score_cor = Biobase::fData(eset_srt)[, c("score_cor", "pval_cor", "insig", "leading_edge"), drop = FALSE],
    sig_score = Biobase::pData(eset_srt)[, "sig_score", drop = FALSE],
    ks = ks_out
  )

  full_heatmap <- NULL
  sig_heatmap <- NULL
  if (isTRUE(make_heatmap_all)) {
    full_heatmap <- spc_heatmap_all(
      eset = eset,
      spc_out = spc_out,
      col_ha = col_ha,
      name = name,
      ...
    )
  }
  if (isTRUE(make_heatmap_sig)) {
    sig_heatmap <- spc_heatmap_sig(
      eset = eset,
      spc_out = spc_out,
      col_ha = col_ha,
      min_sigsize = min_sigsize,
      name = name,
      ...
    )
  }

  return(list(
    score_cor = spc_out$score_cor,
    sig_score = spc_out$sig_score,
    heatmap_all_genes = full_heatmap,
    heatmap_sig_genes = sig_heatmap,
    ks = ks_out
  ))
}
