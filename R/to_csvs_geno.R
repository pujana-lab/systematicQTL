#' Combine markers map and  genotypes dataframe
#' 
#' @param markers_map Dataset with markers
#' @param genotypes Dataset with phenotypes
#' 
#' @return DataFrame for genotypes mapped in CSVS format
#' 
#' 
#' @details 
#' 
#' Markers map file should contain 3 columns:
#' 
#'  - marker: Genomic marker Name (snp id for example)
#'  - chr: Chromosome
#'  - cm: Distance in centimorgans
#'  
#' Genotypes file must be a data frame containing at least a numeric "id" column and the markers in same order than markers map.
#' 
#' @export
#' 
to_csvs_geno <- function(markers_map, genotypes){
  
  markers_map_t = t(markers_map[,-1]) 
  colnames(markers_map_t) = markers_map[,1]
  markers_map_t = cbind(id = c('',''), markers_map_t)
  markers_map_t = data.frame(markers_map_t, stringsAsFactors = FALSE)
  csvs = rbind(markers_map_t, genotypes)  
  return(csvs)
  
}

