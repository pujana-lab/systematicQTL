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


test_that("Should intialize QTL wrapper as S4 class", {

  #Arrange
  qtl = load_cross('./random.genotypes.csvs', './random.pheno.csv')

  #Act
  qtl_wrapper = QtlWrapper(qtl = qtl)

  #Assert
  expect_equal(class(qtl_wrapper@qtl), 'cross2')
  expect_equal(class(qtl_wrapper@phenotypes), c('character'))
  expect_equal(qtl_wrapper@phenotypes, c('signature1', 'signature2', 'signature3'))
  expect_equal(class(qtl_wrapper@map), 'list')
  expect_equal(class(qtl_wrapper@pr), c('calc_genoprob', 'list'))
  expect_equal(class(qtl_wrapper@genescan), c('scan1', 'matrix'))
  expect_equal(class(qtl_wrapper@Xcovar), c('matrix'))


})


test_that("Should set significance_lod manually", {

  #Arrange
  qtl = load_cross('./random.genotypes.csvs', './random.pheno.csv')
  qtl_wrapper = QtlWrapper(qtl = qtl)

  #Act
  qtl_wrapper = set_significance_lod(qtl_wrapper, 4)

  #Assert
  expect_equal(class(qtl_wrapper@pvalues),'matrix')
  expect_equal(as.numeric(qtl_wrapper@pvalues[1,]), c(4,4,4))

})


test_that("Should set significance_lod via permutations", {

  #Arrange
  set.seed(12345)
  qtl = load_cross('./random.genotypes.csvs', './random.pheno.csv')
  qtl_wrapper = QtlWrapper(qtl = qtl)

  #Act
  qtl_wrapper = set_empirical_significance_lod(qtl_wrapper, 100, 0.05, 1)

  #Assert
  expect_equal(class(qtl_wrapper@pvalues),c('summary.scan1perm', 'matrix'))
  expect_equal(round(as.numeric(qtl_wrapper@pvalues[1,]), digits = 2), c(2.85, 2.77, 2.66))

})


test_that("Should compute peaks when find_peaks method is invoked", {

  #Arrange
  set.seed(12345)
  qtl = load_cross('./random.genotypes.csvs', './random.pheno.csv')
  qtl_wrapper = QtlWrapper(qtl = qtl)
  min_lod = 0.4
  qtl_wrapper = set_significance_lod(qtl_wrapper, min_lod)

  #Act
  qtl_wrapper = find_peaks(qtl_wrapper)

  #Assert
  expect_equal(class(qtl_wrapper@peaks),c('data.frame'))
  expect_true(all(qtl_wrapper@peaks$lod > min_lod))

})

