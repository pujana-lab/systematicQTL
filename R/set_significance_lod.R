

#' Sets significance LOD threshold for alpha < 0.05 manually.
#'
#' @param object An QTL object
#' @param min_lod The minimal threshold value for all signatures
#'
#' @returns qtl
#'
#' @export
setGeneric(
  "set_significance_lod",
  function(object, min_lod) {
    standardGeneric("set_significance_lod")
  }
)

#' Sets significance LOD threshold for alpha < 0.05 manually.
#'
#' @param object An QTL object
#' @param min_lod The minimal threshold value for all signatures
#'
#' @returns qtl
#'
#' @export
setMethod("set_significance_lod", 'QtlWrapper', function(object, min_lod){

  pvalues = t(data.frame("0.05" = rep(min_lod, ncol(object@qtl$pheno))))
  colnames(pvalues) = colnames(object@qtl$pheno)
  object@pvalues = pvalues
  object
})


