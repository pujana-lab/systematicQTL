
#' Sets significance LOD threshold empirically using permutation tests
#'
#' @param object An QTL object
#' @param permutations Integer with the number of permutations to be executed
#' @param alpha Float value with the alpha value to be chosen, i.e. 0.05
#' @param cores Number of cores to use for run this permutation test
#'
#' @returns qtl
#'
#' @export
setGeneric(
  "set_empirical_significance_lod",
  function(object, permutations, alpha, cores) {
    standardGeneric("set_empirical_significance_lod")
  }
)



#' Sets significance LOD threshold empirically using permutation tests
#'
#' @param object An QTL object
#' @param permutations Integer with the number of permutations to be executed
#' @param alpha Float value with the alpha value to be chosen, i.e. 0.05
#' @param cores Number of cores to use for run this permutation test
#'
#' @returns qtl
#'
#' @export
setMethod("set_empirical_significance_lod", 'QtlWrapper', function(object, permutations, alpha, cores){

    operm <- qtl2::scan1perm(
      object@pr,
      object@qtl$pheno,
      Xcovar=object@Xcovar,
      n_perm=permutations,
      cores=cores,
      addcovar = object@qtl$covar
    )

    summary(operm, alpha = c(0.2, 0.10, 0.05))

    object@pvalues = summary(operm, alpha = alpha)
    object
})
