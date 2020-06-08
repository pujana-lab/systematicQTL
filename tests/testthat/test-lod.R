source('./load_cross.R')

test_that("Should return LOD scores data frame without p.values when invoked before setting pvalue", {

  #Arrange
  set.seed(12345)
  qtl = load_cross('./random.genotypes.csvs', './random.pheno.csv')
  qtl_wrapper = QtlWrapper(qtl = qtl)

  #Act
  lods = lod(qtl_wrapper)

  qtl_genotypes = as.character(
    do.call('c', lapply(qtl_wrapper@qtl$geno, function(gene){
        return(colnames(gene))
      })
    )
  )

  #Assert
  expect_equal(sort(as.character(lods$snp)), sort(qtl_genotypes))
  expect_equal(colnames(lods), c('snp', 'chr', 'cm',  'signature1', 'signature2', 'signature3'))
})

test_that("Should return LOD scores data frame with p.values when invoked after setting pvalue", {

  #Arrange
  set.seed(12345)
  qtl = load_cross('./random.genotypes.csvs', './random.pheno.csv')
  qtl_wrapper = QtlWrapper(qtl = qtl)

  #Act
  qtl_wrapper = set_significance_lod(qtl_wrapper, 1)
  lods = lod(qtl_wrapper)

  qtl_genotypes = as.character(
    do.call('c', lapply(qtl_wrapper@qtl$geno, function(gene){
        return(colnames(gene))
      })
    )
  )

  #Assert
  expect_equal(sort(as.character(lods$snp)), sort(qtl_genotypes))

  expect_equal(colnames(lods)[1:6],  c('snp', 'chr', 'cm', 'signature1', 'signature2', 'signature3'))
  expect_equal(colnames(lods)[7:8],  c('signature1.threshold', 'signature1.signif'))
  expect_equal(colnames(lods)[9:10],  c('signature2.threshold', 'signature2.signif'))
  expect_equal(colnames(lods)[11:12], c('signature3.threshold', 'signature3.signif'))
})


