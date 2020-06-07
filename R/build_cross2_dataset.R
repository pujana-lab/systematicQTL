#' build_cross2_dataset
#'
#' Generates cross2 dataset from basic geno and phenotypes name. Geno file should combine genenotype values and mapping (see to_csvs_geno) and phenotypes table should be a matrix-type object (users in rows, phenotypes in columns)
#'
#' @param geno_filename Path to csvs geno file
#' @param pheno_filename Path to csv pheno file
#' @param covar_column_names List of "phenotype" columns related to adjusting covariates.
#' @param pheno_column_names List of "phenotype" columns realted to main variables to study.
#' @param genotypes List of allowed genotype values of genotypes file
#' @param alleles  List of homozygote alleles included, should be a subset of genotypes variable
#' @param na.strings List of strings associated to NA values
#' @param comment.char List or string of allowed comment chars
#' @param ... Additional arguments passed to read.cross function
#'
build_cross2_dataset <- function(
  geno_filename,
  pheno_filename,
  covar_column_names,
  pheno_column_names,
  genotypes = c('A', 'H','B'),
  alleles  = c('A', 'B'),
  na.strings = c('-1', 'discrepancy', 'NA'),
  comment.char = '#',
  ...
) {

   dataset = qtl::read.cross(
    'csvs',
    genfile = geno_filename,
    phefile = pheno_filename,
    genotypes = genotypes,
    alleles = alleles,
    na.strings = na.strings,
    comment.char = comment.char,
    convertXdata = FALSE,
    ...
  )

  dataset = qtl::jittermap(dataset)

  dataset = qtl2::convert2cross2(dataset)

  if(is.null(dataset$covar)){
    dataset$covar = data.frame(row.names = rownames(dataset$pheno))
  }

  for(covariate in covar_column_names){
    dataset$covar[,covariate] = dataset$pheno[,covariate]
  }

  dataset$pheno = dataset$pheno[,colnames(dataset$pheno) %in% pheno_column_names]

  return(dataset)

}
