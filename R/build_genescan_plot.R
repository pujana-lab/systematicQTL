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
