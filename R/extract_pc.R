#' Extract a principal component score for a signature
#'
#' @param eset an \code{ExpressionSet}
#' @param sig character vector of feature names in \code{eset}
#' @param i index of the principal component to extract
#'
#' @return named numeric vector of PC scores (one per sample)
#' @importFrom methods is
#' @importFrom Biobase ExpressionSet exprs featureNames sampleNames
#'
#' @export
extract_pc <- function(eset, sig, i = 1) {
  stopifnot(methods::is(eset, "ExpressionSet"))
  stopifnot(is.character(sig), length(sig) > 0)
  stopifnot(length(i) == 1, is.numeric(i), i == as.integer(i), i >= 1)
  stopifnot(all(sig %in% Biobase::featureNames(eset)))

  dat <- Biobase::exprs(eset)[sig, , drop = FALSE]
  pca <- stats::prcomp(t(dat), center = TRUE, scale. = FALSE)
  if (i > ncol(pca$x)) {
    stop("requested PC index exceeds available components")
  }
  pc <- pca$x[, i]
  names(pc) <- Biobase::sampleNames(eset)
  return(pc)
}
