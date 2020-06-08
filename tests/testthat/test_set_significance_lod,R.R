source('./load_cross.R')
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

