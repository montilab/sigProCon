## Internal helpers shared across R/ files. Not exported.

## Coerce a SummarizedExperiment to an ExpressionSet; ExpressionSet input
## passes through unchanged. colData()/rowData() return an S4 DFrame, which
## Biobase::AnnotatedDataFrame() has no method for, hence the as.data.frame().
#' @keywords internal
.as_expressionset <- function(eset) {
  if (methods::is(eset, "SummarizedExperiment")) {
    eset <- Biobase::ExpressionSet(
      assayData = SummarizedExperiment::assay(eset),
      phenoData = Biobase::AnnotatedDataFrame(as.data.frame(SummarizedExperiment::colData(eset))),
      featureData = Biobase::AnnotatedDataFrame(as.data.frame(SummarizedExperiment::rowData(eset)))
    )
  }
  eset
}

## Warn when a column this package is about to write already exists in the
## caller's pData/fData, since its prior values are about to be overwritten.
#' @keywords internal
.warn_if_overwriting <- function(existing_names, col) {
  if (col %in% existing_names) {
    warning(
      "column '", col, "' already exists in the input eset and will be overwritten",
      call. = FALSE
    )
  }
  invisible(NULL)
}
