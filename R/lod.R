
#' Retrieves real values LOD scores from a QTL wrapper object. If P.values have been defined (manually or empirically) also
#' returns significance stores
#'
#' @param object An QTL object
#'
#' @returns Data frame with genotype LOD scores
#'
#' @export
setGeneric(
  "lod",
  function(object) {
    standardGeneric("lod")
  }
)


#' Retrieves real values LOD scores from a QTL wrapper object. If P.values have been defined (manually or empirically) also
#' returns significance stores
#'
#' @param object An QTL object
#'
#' @returns Data frame with genotype LOD scores
#'
#' @export
setMethod("lod", 'QtlWrapper', function(object){

  snp_lods = as.data.frame(
    object@genescan[!grepl('\\.loc', rownames(object@genescan)),]
  )

  map_markers = as.data.frame(unlist(object@qtl$gmap))
  map_markers = data.frame(
    snp = gsub('^[[:alnum:]]+\\.', '', rownames(map_markers)),
    chr = gsub('\\..*', '', rownames(map_markers)),
    cm  = map_markers[,1]
  )

  snp_lods = merge(
    map_markers,
    snp_lods,
    by.x = 'snp',
    by.y = 0
  )

  if(!is.null(object@pvalues)){
    for(pvalueName in colnames(object@pvalues)){
      snp_lods[,sprintf('%s.threshold', pvalueName)] = object@pvalues[,pvalueName]
      snp_lods[,sprintf('%s.signif', pvalueName)] = snp_lods[,pvalueName] > object@pvalues[,pvalueName]
    }
  }
  return(snp_lods)
})


