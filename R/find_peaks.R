#' Finds the peaks passing lod values.
#'
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
