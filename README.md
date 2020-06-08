# Systematic QTL

This package contains the software tools developed for running the experiments which support the QTL-related results explained at "Immune Cell Associations with Cancer Risk", paper pending of publication. 

This is a systematic approach for [QTL pipeline developed and maintained by Karl Broman](https://github.com/rqtl/qtl2/).


## Installation
  
To install it, the easiest is to use the R package devtools and its function install_github. To do so, open an R session and enter

```
install.packages(c("devtools","curl")) ##Installs devtools 
library(devtools)
install_github("pujana-lab/systematicQTL",ref="master")
```

## Basic running

### File description

This pipeline requires to be run 3 files: 
 - Phenotypes data table in matrix format, where phenotype values mare in columns and users in rows. This file must be comma separated and first column (id column)  should be the user identifier.
 - A file with this tree columns:
    - snp: The snp identifier
    - chr: The chromosome number (1-22, X, Y) where the SNP is located.
    - cm:  SNP distance (in centimorgans) to begining of the chromosome
 - A genotype data table file, where first column (id) is for user id, and other values are genotype values.



### Object definition

Due R/QTL2 file requires filenames to be build, is necessary construct this object in two steps: 
 - Define the genotype CSVS file from genotype and mapping files.
 - Call wrapper constructor setting also main phenotype column names and genotype ones

For example.

```
geno_filename = '../tests/testthat/random.genotypes.csvs'
pheno_filename  = '../tests/testthat/random.pheno.csv'

chosen_alleles      = c('A','B')
chosen_covars       = c('pheno_var1', 'pheno_var2', 'pheno_var3')
chosen_phenotypes   = c('signature1', 'signature2', 'signature3')

cross2 = systematicQTL::build_cross2_dataset(
  geno_filename = geno_filename,
  pheno_filename = pheno_filename,
  covar_column_names = chosen_covars,
  pheno_column_names = chosen_phenotypes,
  alleles = chosen_alleles
)

qtl_wrapper = systematicQTL::QtlWrapper(qtl = cross2)
```
 
### LOD thresholds definition


There are to ways to define thresholds to considere LOD scores significant, the thresholds can be defined manually or using a [permutation test](https://kbroman.org/qtl2/assets/vignettes/user_guide.html#performing_a_permutation_test). 

```
qtl_wrapper = systematicQTL::set_significance_lod(qtl_wrapper,0.5)

qtl_wrapper = systematicQTL::set_empirical_significance_lod(qtl_wrapper, 100, 0.05, 1)
```


### Peaks identification

Next step is identify the associated peaks using [find_peaks](https://kbroman.org/qtl2/assets/vignettes/user_guide.html#finding_lod_peaks) considering previously computed threshold. Also,  helping method build_genescan_plot allows describe the different peak plots.

```
qtl_wrapper = systematicQTL::find_peaks(qtl_wrapper)

for(phenotype in qtl_wrapper@phenotypes){
  systematicQTL::build_genescan_plot(qtl_wrapper, phenotype, sprintf('Pheno: %s', phenotype))
}

```
 
### SNP assocaition and LOD score tables

`Lod` method returns the LOD tables for all snps included annotating possible significant correlations with profided signatures. To describe them method `build_pdx_plot` can be executed, as in below example.

```

lods = systematicQTL::lod(qtl_wrapper)

for( i in 1:nrow(lods)){
  this_lod = lods[i,]
  systematicQTL::build_pdx_plot(
    qtl_wrapper, chr = this_lod$chr, cm = this_lod$cm, 'signature1', 
    sprintf('%s-signature1', this_lod$snp),
    sprintf('LOD score: %.3f, (tresh. %.3f)', this_lod$signature1, this_lod$signature1.threshold)
  )
}



```
 
