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
    file = './random.pheno.csv',
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
  write_csvs_geno(csvs, './random.genotypes.csvs')
}