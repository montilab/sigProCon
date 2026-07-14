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
