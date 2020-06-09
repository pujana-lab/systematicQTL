# Identification of immune/stromal cell quantitative trait loci linked to cancer risk

This package contains the pipeline and tools developed for the study entitled "Immune Cell Associations with Cancer Risk". This is pipeline was developed and is maintained by Luis Palomero and Roderic Espin (MA Pujana’s lab, Catalan Institute of Oncology, IDIBELL).

The pipeline includes the R package ConsensusTME (Jiménez-Sánchez et al., 2019; https://github.com/cansysbio/ConsensusTME), ssGSEA in Gene Set Variation Analysis (GSVA) (Hänzelmann et al., 2013; 10.18129/B9.bioc.GSVA), R/qtl2 (Broman et al., 2019; https://github.com/rqtl/qtl2), and bestNormalize (https://github.com/petersonR/bestNormalize).



## Installation
  
To install it use the R package devtools and its function install_github. Open an R session and enter the following commands:

```
install.packages(c("devtools","curl")) ##Installs devtools 
library(devtools)
install_github("pujana-lab/systematicQTL",ref="master")
```

## Basic running

### File description

This pipeline requires to run 3 files: 
 - Phenotypes data table in matrix format, where phenotype values are in columns and cases in rows. This file must be comma separated and first column (ID column) should be case identifier.
 - A file with these three columns:
    - snp: snp identifier
    - chr: chromosome number (1-22, X, Y) where SNP maps.
    - cm:  SNP distance (cM) to chromosome start.
 - A genotype data table, where first column (ID) is for cases and other values are genotypes.


### Object definition

Given that R/QTL2 requires filenames to be build, it is necessary to generate this object in two steps: 
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


There are two ways to define LOD thresholds, manually or using permutations.

```
qtl_wrapper = systematicQTL::set_significance_lod(qtl_wrapper,0.5)

qtl_wrapper = systematicQTL::set_empirical_significance_lod(qtl_wrapper, 100, 0.05, 1)
```


### Peaks identification

Next step is to identify the associated peaks using find_peaks considering previously computed thresholds. Also, build_genescan_plot allows to obtain peak plots.

```
qtl_wrapper = systematicQTL::find_peaks(qtl_wrapper)

for(phenotype in qtl_wrapper@phenotypes){
  systematicQTL::build_genescan_plot(qtl_wrapper, phenotype, sprintf('Pheno: %s', phenotype))
}

```
 
### SNP association and LOD score tables

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
 
