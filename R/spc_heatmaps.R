#' Generate All-Genes SPC Heatmap
#'
#' Build the full (all-genes) heatmap from a
#' \code{signature_projection_contributors()} result.
#'
#' @param eset an expression set
#' @param spc_out output list from \code{signature_projection_contributors()}
#' @param col_ha a ComplexHeatmap::heatmapAnnotation object with columns' (i.e., samples') annotation
#' @param name name for the heatmap
#' @param row_title row title passed to \code{ComplexHeatmap::Heatmap}
#' @param show_row_names logical, whether to show row names
#' @param show_column_names logical, whether to show column names
#' @param row_names_side side for row names (\code{"left"} or \code{"right"})
#' @param ... additional parameters to pass to ComplexHeatmap::Heatmap
#'
#' @return a ComplexHeatmap object
#'
#' @import Biobase
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation anno_barplot
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom methods is
#' @importFrom grid gpar
#'
#' @export
spc_heatmap_all <- function(
    eset,
    spc_out,
    col_ha = NULL,
    name = "expression",
    ## ComplexHeatmap::Heatmap arguments
    row_title = "Genes",
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_names_side = "left",
    fontsize = 8,
    ...
) {
  prep <- .spc_prepare_heatmap_data(eset = eset, spc_out = spc_out, col_ha = col_ha)

  row_ha <- ComplexHeatmap::rowAnnotation(
    genes = prep$score_cor$insig,
    leadedge = ifelse(rownames(prep$score_cor) %in% prep$ks_hits, "yes", "no"),
    correlation = ComplexHeatmap::anno_barplot(prep$score_cor$score_cor),
    col = list(
      genes = c("background" = "brown", "signature" = "lightgreen"),
      leadedge = c(yes = "black", no = "white")
    ),
    show_annotation_name = FALSE
  )
  suppressMessages(ComplexHeatmap::Heatmap(
    matrix = t(scale(t(Biobase::exprs(prep$eset_srt)))),
    top_annotation = prep$col_ha_srt,
    name = name,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    row_title = row_title,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_names_side = row_names_side,
    column_names_gp = grid::gpar(fontsize = fontsize),
    ...
  )) + row_ha
}

#' Generate Signature-Only SPC Heatmap
#'
#' Build the signature-only heatmap from a
#' \code{signature_projection_contributors()} result.
#'
#' @param eset an expression set
#' @param spc_out output list from \code{signature_projection_contributors()}
#' @param col_ha a ComplexHeatmap::heatmapAnnotation object with columns' (i.e., samples') annotation
#' @param min_sigsize minimum number of signature genes required to plot signature-only heatmap
#' @param name name for the heatmap
#' @param row_title row title passed to \code{ComplexHeatmap::Heatmap}
#' @param show_row_names logical, whether to show row names
#' @param show_column_names logical, whether to show column names
#' @param row_names_side side for row names (\code{"left"} or \code{"right"})
#' @param ... additional parameters to pass to ComplexHeatmap::Heatmap
#'
#' @return a ComplexHeatmap object
#'
#' @import Biobase
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation anno_barplot
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom methods is
#' @importFrom grid gpar
#'
#' @export
spc_heatmap_sig <- function(
    eset,
    spc_out,
    col_ha = NULL,
    min_sigsize = 3,
    name = "expression",
    ## ComplexHeatmap::Heatmap arguments
    row_title = "Signature Genes",
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_side = "left",
    fontsize = 8,
    ...
) {
  prep <- .spc_prepare_heatmap_data(eset = eset, spc_out = spc_out, col_ha = col_ha)
  sig_genes <- rownames(prep$score_cor)[prep$score_cor$insig == "signature"]
  eset_flt <- prep$eset_srt[sig_genes, ]
  score_cor_sig <- prep$score_cor[sig_genes, "score_cor", drop = TRUE]
  stopifnot(nrow(eset_flt) > min_sigsize)

  suppressMessages(ComplexHeatmap::Heatmap(
    matrix = t(scale(t(Biobase::exprs(eset_flt)))),
    top_annotation = prep$col_ha_srt,
    name = name,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    row_title = row_title,
    row_title_gp = grid::gpar(fontsize = fontsize*1.5, fontface = "bold"),
    row_names_gp = grid::gpar(fontsize = fontsize),
    column_names_gp = grid::gpar(fontsize = fontsize),
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_names_side = row_names_side,
    ...
  )) + ComplexHeatmap::rowAnnotation(
    correlation = ComplexHeatmap::anno_barplot(score_cor_sig)
  )
}

.spc_prepare_heatmap_data <- function(eset, spc_out, col_ha = NULL) {
  stopifnot(methods::is(spc_out, "list"))
  stopifnot(all(c("score_cor", "sig_score") %in% names(spc_out)))
  stopifnot(methods::is(eset, "SummarizedExperiment") || methods::is(eset, "ExpressionSet"))

  if (methods::is(eset, "SummarizedExperiment")) {
    eset <- Biobase::ExpressionSet(
      assayData = SummarizedExperiment::assay(eset),
      phenoData = Biobase::AnnotatedDataFrame(SummarizedExperiment::colData(eset)),
      featureData = Biobase::AnnotatedDataFrame(SummarizedExperiment::rowData(eset))
    )
  }

  score_cor <- spc_out$score_cor
  sig_score <- spc_out$sig_score
  stopifnot(methods::is(score_cor, "data.frame"))
  stopifnot(methods::is(sig_score, "data.frame"))
  stopifnot(all(c("score_cor", "insig") %in% colnames(score_cor)))
  stopifnot("sig_score" %in% colnames(sig_score))
  stopifnot(all(rownames(score_cor) %in% Biobase::featureNames(eset)))
  stopifnot(all(rownames(sig_score) %in% Biobase::sampleNames(eset)))

  row_ord <- rownames(score_cor)
  col_ord <- rownames(sig_score)
  eset_srt <- eset[row_ord, col_ord]

  ks_hits <- character(0)
  if (!is.null(spc_out$ks) && !is.null(spc_out$ks$hits)) {
    if (!all(is.na(spc_out$ks$hits))) {
      ks_hits <- spc_out$ks$hits
    }
  }

  if (is.null(col_ha)) {
    col_ha_srt <- ComplexHeatmap::HeatmapAnnotation(
      sig_score = ComplexHeatmap::anno_barplot(sig_score[col_ord, "sig_score"])
    )
  } else {
    stopifnot(isTRUE(all(rownames(col_ha) %in% Biobase::sampleNames(eset))))
    col_ha_srt <- c(
      ComplexHeatmap::HeatmapAnnotation(
        sig_score = ComplexHeatmap::anno_barplot(sig_score[col_ord, "sig_score"])
      ),
      col_ha[match(col_ord, Biobase::sampleNames(eset)), ]
    )
  }

  list(
    eset_srt = eset_srt,
    score_cor = score_cor[row_ord, , drop = FALSE],
    col_ha_srt = col_ha_srt,
    ks_hits = ks_hits
  )
}
