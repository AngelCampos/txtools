
<!-- README.md is generated from README.Rmd. Please edit that file -->

# txtools

<!-- badges: start -->

<!-- badges: end -->

The goal of txtools is to process paired-end reads data from a genomic space 
into their transcriptomic versions. In doing so we need to allocate both reads
into their corresponding gene model(s) resulting in a unified read. This 
representation allows us to manipulate with greater confidence and ease 
paired-end reads which may be skipping exons and/or transversing long intronic 
regions.

## Installation

Install txtools dependencies with the following code

``` r
# CRAN packages
install.packages(c("magrittr", "stringr", "devtools"))
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "plyranges", "Rsamtools", "GenomicAlignments"))
```

Install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AngelCampos/testRepo")
BiocManager::install("AngelCampos/txtools")
```

## Demo

This is a basic example which shows you how to solve a common problem:

``` r
# library(txtools)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
