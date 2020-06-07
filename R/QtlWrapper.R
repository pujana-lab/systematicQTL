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


#' Finds the peaks passing lod values.
#' This method should be called after defining LOD values manually or empirically
#'
#' @param object A QTL object
#'
#' @returns qtl
#'
#' @export
setGeneric(
  "find_peaks",
  function(object){
    standardGeneric("find_peaks")
  }
)


#' Finds the peaks passing lod values.
#' This method should be called after defining LOD values manually or empirically
#'
#' @param object A QTL object
#'
#' @returns qtl
#'
#' @export
setMethod("find_peaks", "QtlWrapper", function(object){
  object@peaks = qtl2::find_peaks(
    object@genescan,
    object@map,
    threshold=object@pvalues,
    expand2markers = TRUE,
    drop = 0.5
  )
  object
})


#' Draw genescan peaks plot
#'
#' @param object A QTL object
#' @param phenotype Phenotype label to be displayed.
#' @param title Main title value
#'
#' @returns plot
#'
#' @export
setGeneric(
  "build_genescan_plot",
  function(object, phenotype, title){
    standardGeneric("build_genescan_plot")
  }
)


#' Draw genescan peaks plot
#'
#' @param object A QTL object
#' @param phenotype Phenotype label to be displayed.
#' @param title Main title value
#'
#' @returns plot
#'
#' @export
setMethod("build_genescan_plot", "QtlWrapper", function(object, phenotype, title){

  graphics::plot(
    object@genescan,
    object@map,
    lodcolumn=phenotype, col="slateblue",
    main = title

  )
  graphics::legend(
    "topright",
    lwd=2,
    col=c("slateblue"),
    phenotype,
    bg="gray90",
    cex = .8,
  )
  graphics::abline(h = object@pvalues[,phenotype], col = 'darkred', lty = 2)
})

