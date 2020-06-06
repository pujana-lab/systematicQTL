#' Helper function to export csvs_geno dataset in correct format expected as input for buildCross2Dataset
#' 
#' @param csvs_geno Dataframe (output of build_genos_df)
#' @param filename Where the data would be
#' 
#' @return NULL
#' 
#' @export
#' 
write_csvs_geno <- function(csvs_geno, filename){
  utils::write.table(
    csvs_geno,
    sep = ',',
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    file = filename,
    na = ''
  )
}