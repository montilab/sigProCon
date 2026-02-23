#####################################################
## EXTRACT EIGENGENES
#####################################################
## given an ExpressionSet and a list of "signatures" (sets of features), 
## ..extract the eigengene (i.e., 1st PCA) from each module and 
## ..return as an ExpressionSet
##
#' Extract module eigengenes from an expression set
#'
#' @param eset an \code{ExpressionSet}
#' @param signatures named list of character vectors of features
#' @param mod_fdata optional feature-level metadata for modules
#' @param key if \code{NULL} match signatures to \code{featureNames(eset)},
#'   otherwise use \code{fData(eset)[, key]}
#' @param min_size minimum number of overlapping features required per signature
#' @param method eigengene extraction method; use \code{"single"} to run modules
#'   independently and \code{"multiple"} for one joint call
#' @param name_cleanup logical; remove leading \code{"ME"} from module names
#'
#' @return an \code{ExpressionSet} with eigengenes as features and input samples as
#'   columns
#'
#' @export
extract_eigengenes <- function(
    eset,
    signatures,       # p-sized list of character vectors
    mod_fdata = NULL, # pxq data frame of signatures annotations
    key = NULL,       # if NULL use featureNames, otherwise use fData(eset)[,key]
    min_size = 3,
    method = c("multiple", "single"),
    name_cleanup = TRUE) 
{
  ## checks
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("package 'purrr' is required for extract_eigengenes()")
  }
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("package 'WGCNA' is required for extract_eigengenes()")
  }
  method <- match.arg(method)
  if (class(eset) != "ExpressionSet") {
    stop("ExpressionSet expected:", class(eset))
  }
  if (is.null(key) && length(intersect(Biobase::featureNames(eset), unique(unlist(signatures)))) < 1) {
    stop("no overlap between eset features and signatures")
  }
  if (!is.null(key) && !(key %in% colnames(Biobase::fData(eset)))) {
    stop("unrecognized fData key:", key)
  }
  if (!is.null(mod_fdata) && !all(names(signatures) %in% rownames(mod_fdata))) {
    stop("!all(names(signatures) %in% rownames(mod_fdata))")
  }
  if (method == "multiple" && length(signatures) > 1) {
    tmp <- utils::combn(signatures, 2, function(mod) length(intersect(mod[[1]], mod[[2]])), simplify = FALSE)
    if (any(tmp > 0)) stop("method=='multiple' can't be used with overlapping signatures")
  }
  ## end checks

  ## for easier handling, create new fdata column
  Biobase::fData(eset)$key <- if (is.null(key)) Biobase::featureNames(eset) else Biobase::fData(eset)[, key]

  ## remove signatures w/ not enough features
  signatures1 <- signatures[vapply(signatures, function(Z) {
    length(intersect(Z, Biobase::fData(eset)$key))
  }, numeric(1)) >= min_size]
  if (length(signatures1) < 1) stop("no module left")

  ## reduce eset to features in signatures
  eset1 <- eset[Biobase::fData(eset)$key %in% unique(unlist(signatures1)), ]

  ## calculate eigengene for each module (map2_dfc iterates over pairs of input args)
  if (method == "single") {
    eigenDat <- purrr::map2_dfc(signatures1, names(signatures1), function(module, modname) {
      datI <- t(Biobase::exprs(eset1)[Biobase::fData(eset1)$key %in% module, ])
      WGCNA::moduleEigengenes(expr = datI, colors = rep(modname, ncol(datI)))$eigengenes
    })
  }
  ## this is only guaranteed to work with non-overlapping signatures
  else if (method == "multiple") {
    ## probably there's a more elegant way to generate the COL vector
    COL <- factor(rep(NA, nrow(eset1)), levels = sort(names(signatures1)))
    for (i in seq(signatures1)) COL[Biobase::fData(eset1)$key %in% signatures1[[i]]] <- names(signatures1)[i]
    eigenDat <- WGCNA::moduleEigengenes(expr = t(Biobase::exprs(eset1)), colors = COL, excludeGrey = TRUE)$eigengenes
  } else {
    stop("unrecognized method: ", method)
  }
  ## cleanup module names
  if (name_cleanup) {
    colnames(eigenDat) <- gsub("\\bME", "", colnames(eigenDat))
  }
  ## repackage in an ExpressionSet
  fdata <- {
    if (is.null(mod_fdata)) {
      data.frame(module = names(signatures1), row.names = colnames(eigenDat))
    } else {
      mod_fdata[colnames(eigenDat), , drop = FALSE]
    }
  }
  eigenSet <- Biobase::ExpressionSet(
    assay = as.matrix(t(eigenDat)),
    phenoData = Biobase::AnnotatedDataFrame(Biobase::pData(eset1)),
    featureData = Biobase::AnnotatedDataFrame(fdata)
  )
  return(eigenSet)
}
