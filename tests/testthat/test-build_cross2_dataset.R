

context("Build_cross2_dataset")

{
  set.seed(12345)
  n_samples = 100
  pheno_data = data.frame(
    id = 1:n_samples,
    pheno_var1 = sample(0:3, n_samples, replace = TRUE),
    pheno_var2 = sample(0:4, n_samples, replace = TRUE),
    pheno_var3 = sample(0:4, n_samples, replace = TRUE),
    sex =        sample(0:1, n_samples, replace = TRUE),
    signature1 = rnorm(n_samples, 0, 3),
    signature2 = rnorm(n_samples, 10, 2),
    signature3 = rexp(n_samples, 1)
  )

  write.table(
    pheno_data,
    sep = ',',
    file = './pheno.csv',
    row.names = FALSE,
    col.names = TRUE
  )


  n_markers = 15

  map_data = data.frame(
    snp = sprintf('rs0000%s', 1:n_markers),
    chr = sample(c(1:22, 'X'), n_markers, replace = TRUE),
    cm = runif(n_markers, 1, 100)
  )


  genos = data.frame(
    id = 1:n_samples
  )

  for(marker in map_data$snp){
    marker_chr = map_data[map_data$snp == marker, 'chr']
    if(marker_chr == 'X'){
      genos[,marker] = sample(c('A','H'), n_samples,  replace = TRUE)
    } else {
      genos[,marker] = sample(c('A','H','B'), n_samples,  replace = TRUE)
    }
  }

  csvs = to_csvs_geno(map_data, genos)
  write_csvs_geno(csvs, './genos.csvs')
}

test_that("multiplication works", {

  #Arrange
  geno_filename = './genos.csvs'
  pheno_filename  = './pheno.csv'

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
