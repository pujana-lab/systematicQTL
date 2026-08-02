# systematicQTL

**Identification of immune/stromal cell quantitative trait loci linked to cancer risk**

This package contains the pipeline and tools developed for the study:

**Palomero L, Galván-Femenía I, de Cid R, Espín R, Barnes DR, CIMBA, et al. Antoniou AC, Lázaro C, Pujana MA. Immune Cell Associations with Cancer Risk.** *iScience*. 2020;23(7):101296. [DOI: 10.1016/j.isci.2020.101296](https://doi.org/10.1016/j.isci.2020.101296) · [PMID: 32622267](https://pubmed.ncbi.nlm.nih.gov/32622267)

This pipeline was developed and is maintained by Luis Palomero and Roderic Espín (MA Pujana's lab, Catalan Institute of Oncology, IDIBELL).

> **Note:** This package performs QTL mapping (LOD scoring, peak identification). For the upstream data pre-processing and regression step (case filtering, normalization, covariate selection), see [systematicBNR](https://github.com/pujana-lab/systematicBNR).

The pipeline includes the R package ConsensusTME (Jiménez-Sánchez et al., 2019; [github.com/cansysbio/ConsensusTME](https://github.com/cansysbio/ConsensusTME)), ssGSEA in Gene Set Variation Analysis (GSVA) (Hänzelmann et al., 2013; `10.18129/B9.bioc.GSVA`), R/qtl2 (Broman et al., 2019; [github.com/rqtl/qtl2](https://github.com/rqtl/qtl2)), and bestNormalize ([github.com/petersonR/bestNormalize](https://github.com/petersonR/bestNormalize)).

## Installation

To install it use the R package devtools and its function install_github. Open an R session and enter the following commands:

```r
install.packages(c("devtools","curl")) ##Installs devtools
library(devtools)
install_github("pujana-lab/systematicQTL",ref="master")
```

## Basic running

### File description

This pipeline requires 3 files:

- Phenotypes data table in matrix format, where phenotype values are in columns and cases in rows. This file must be comma separated and the first column (ID column) should be the case identifier.
- A file with these three columns:
  - `snp`: SNP identifier
  - `chr`: chromosome number (1–22, X, Y) where the SNP maps
  - `cm`: SNP distance (cM) to chromosome start
- A genotype data table, where the first column (ID) is for cases and other values are genotypes.

### Object definition

Given that R/QTL2 requires filenames to be built, it is necessary to generate this object in two steps:

- Define the genotype CSVS file from genotype and mapping files.
- Call the wrapper constructor, setting the main phenotype column names and genotype ones.

Example:

```r
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

There are two ways to define LOD thresholds: manually, or using permutations.

```r
qtl_wrapper = systematicQTL::set_significance_lod(qtl_wrapper,0.5)

qtl_wrapper = systematicQTL::set_empirical_significance_lod(qtl_wrapper, 100, 0.05, 1)
```

### Peaks identification

The next step is to identify the associated peaks using `find_peaks`, considering the previously computed thresholds. `build_genescan_plot` allows you to obtain peak plots.

```r
qtl_wrapper = systematicQTL::find_peaks(qtl_wrapper)

for(phenotype in qtl_wrapper@phenotypes){
  systematicQTL::build_genescan_plot(qtl_wrapper, phenotype, sprintf('Pheno: %s', phenotype))
}
```

### SNP association and LOD score tables

The `lod` method returns the LOD tables for all SNPs, annotating possible significant correlations with the provided signatures. `build_pdx_plot` can be used to describe them, as in the example below.

```r
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

## Citation

If you use this code, please cite the publication above.

## Authors

Luis Palomero and Roderic Espín
