source("./load_cross.R")
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
  expect_equal(colnames(qtl_wrapper@Xcovar), c('sex'))

})


