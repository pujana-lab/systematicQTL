load_cross <- function(geno_filename, pheno_filename){

  geno_filename = './random.genotypes.csvs'
  pheno_filename  = './random.pheno.csv'


  chosen_alleles      = c('A','B')
  chosen_covars       = c('pheno_var1', 'pheno_var2', 'pheno_var3')
  chosen_phenotypes   = c('signature1', 'signature2', 'signature3')

  #Act

  cross2 = build_cross2_dataset(
    geno_filename = geno_filename,
    pheno_filename = pheno_filename,
    covar_column_names = chosen_covars,
    pheno_column_names = chosen_phenotypes,
    alleles = chosen_alleles
  )
  return(cross2)
}
