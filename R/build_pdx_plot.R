#' Draw genescan peaks plot
#'
#' @param object A QTL object
#' @param chr The SNP chromosome
#' @param cm  The SNP distance in centimorgans
#' @param phenotype String with phenotype name to display over genotypes
#' @param title String with title plot
#' @param subtitle String with subtitle plot
#'
#' @returns plot
#'
#' @export
setGeneric(
  "build_pdx_plot",
  function(object, chr, cm, phenotype, title, subtitle){
    standardGeneric("build_pdx_plot")
  }
)

#' Draw genescan peaks plot
#'
#' @param object A QTL object
#' @param chr The SNP chromosome
#' @param cm  The SNP distance in centimorgans
#' @param phenotype String with phenotype name to display over genotypes
#' @param title String with title plot
#' @param subtitle String with subtitle plot
#'
#' @returns plot
#'
#' @export
setMethod("build_pdx_plot", "QtlWrapper", function(object, chr, cm, phenotype, title, subtitle){
  g <- qtl2::maxmarg(object@pr, object@map, chr=chr, pos=cm, return_char=TRUE)
  qtl2::plot_pxg(g, object@qtl$pheno[,phenotype], ylab=phenotype,
     main = title,
     sub = subtitle,
     sort = FALSE,
     SEmult = 2
    )
  }
)
