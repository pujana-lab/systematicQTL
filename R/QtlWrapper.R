#' Qtl Object wrapper class.
#'
#'
#' @slot qtl QTL object
#' @slot phenotypes Phenotype names
#' @slot map Map object generated automatically at initialization
#' @slot pr Calc_genoprob object generated automatically at initialization
#' @slot Xcovar Matrix object with Xprobatilities
#' @slot genescan  A scan1 object generated at initialiation
#' @slot pvalues Data frame with minimal LOD values thresholds by signature that should be considered as significant.
#' @slot peaks Data Frame with computed peaks. Set at "find_peaks" method
#'
#' @export
setClass( "QtlWrapper",
  slots = list(
    qtl = 'ANY',
    phenotypes = 'character',
    map = 'ANY',
    pr = 'ANY',
    Xcovar = 'ANY',
    genescan = 'ANY',
    pvalues = 'ANY',
    peaks = 'ANY'
  ))

#' QtlWrapper constructor
#' @param qtl An QTL object
#'
#' @returns qtl
#'
#' @export
QtlWrapper <- function(qtl) {

  map = qtl2::insert_pseudomarkers(qtl$gmap, step=1)
  pr <- qtl2::calc_genoprob(qtl, map, error_prob=0.002)
  Xcovar <- qtl2::get_x_covar(qtl)
  genescan <- qtl2::scan1(pr, qtl$pheno, Xcovar=Xcovar, addcovar = qtl$covar)
  phenotypes = colnames(qtl$pheno)
  methods::new("QtlWrapper",
               qtl = qtl, map = map, pr = pr,
               Xcovar = Xcovar, genescan = genescan,
               phenotypes = phenotypes)
}



