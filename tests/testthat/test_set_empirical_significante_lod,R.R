source("./load_cross.R")
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
