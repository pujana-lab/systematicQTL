source("./load_cross.R")
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

