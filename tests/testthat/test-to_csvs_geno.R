
build_markers_df = function(n_markers){
  data.frame(
    marker = sprintf('rs0000%s', 1:n_markers),
    chr = sample(c(1:22,'X'), n_markers, replace = TRUE),
    cm = runif(n_markers, 1, 100)
  )

}

build_genos_df = function(markers, n_samples){
  genos = data.frame(
    id = 1:n_samples
  )

  for(marker in markers){
    genos[,marker] = sample(c('A','H','B'), n_samples,  replace = TRUE)
  }
  return(genos)
}

test_that("Build csvs genotype file from base geno dataframe and marker map", {

  #Arrange
  n_markers = 25
  n_samples = 100
  markers_map = build_markers_df(n_markers)
  genos = build_genos_df(markers_map$marker, n_samples)

  #Act
  csvs_geno = to_csvs_geno(markers_map, genos)
  #Assert
  expect_equal(nrow(csvs_geno), nrow(genos)+ ncol(markers_map) - 1)
  expect_equal(csvs_geno$id[1:2], c('', ''))

})
