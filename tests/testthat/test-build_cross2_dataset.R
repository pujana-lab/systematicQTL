

context("Build_cross2_dataset")


test_that("multiplication works", {

  #Arrange
  geno_filename = './random.genotypes.csvs'
  pheno_filename  = './random.pheno.csv'

  geno = as.data.frame(data.table::fread(geno_filename))
  pheno = as.data.frame(data.table::fread(pheno_filename))

  chromosomes         = unique(as.character(geno[1,]))
  chromosomes         = chromosomes[chromosomes != 'NA']
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

  #Assert

  expect_equal(class(cross2), 'cross2')
  expect_equal(length(cross2$geno), length(chromosomes))
  expect_equal(cross2$alleles, chosen_alleles)
  expect_equal(colnames(cross2$covar), c('sex', chosen_covars))
  expect_equal(colnames(cross2$pheno), chosen_phenotypes)

})
