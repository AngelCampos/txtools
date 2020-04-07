
# txtools

<!-- badges: start -->

<!-- badges: end -->

**txtools** is a package that processes GenomicAlignments objects into
their transcriptomic versions.

**txtools** is meant to expand the functionality of the
[**GenomicAlignments**](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
package, as currently it does not support transcriptomic convertion.
This type of convertion is increasingly needed to process and analyze
RNA-seq data in which transcript-structure awareness and close
nucleotide-level inspection is required.

## Installation

Install the development version from [GitHub](https://github.com/) with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("AngelCampos/txtools")
```

## Demo

### Starting point `tx_reads()`

The main input that we want to process are Genomic Alignments from **BAM
files**, into their transcriptomic homologues. To do this we require
gene models provided in the form of **BED-12 files**.

In this basic example we use data provided within **txtools**.

We first load the BAM file and the BED file’s gene models.

``` r
# Load packages
library(txtools)

# This example files are installed along txtools
bamFile <- system.file("extdata", "example_hg19.bam", package = "txtools")
bedFile <- system.file("extdata", "twoUCSCgenes_hg19.bed", package = "txtools")

# Loading files and processing them using the gene models
reads <- tx_load_bam(bamFile, loadSeq = T, verbose = T, yieldSize = 10000)
#> Reading number of records in file 
#> 5061 number of BAM records 
#> Loading BAM file 
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#>  
#> 2530 reads succesfully loaded 
#> Dumped reads due to ambiguous pairs: 0
geneAnnot <- tx_load_bed(bedFile) # plyranges read_bed function
```

Then we use the `tx_reads()` function to process the genomic alignments
into transcriptomic versions. txtools uses will assign mappings to their
corresponding genes only if they are overlapping their genomic regions
AND if they are consistent with the exon structure of the gene model,
perfectly distinguishing isoforms.

Currently that means that txtools is designed for and requires RNA-seq
libraries that are **strand-specific**.

**To accelerate processing the multi-core function** `tx_reads_mc()`
**is available for UNIX systems**

To control for spurious mappings we can filter for too long mappings
with the `tx_filter_max_width()` function. In this example we set the
threshold to 500, removing mappings longer than 500 nucleotides at the
transcript level.

``` r
txReads <- tx_reads(reads, geneAnnot, withSeq = T, verbose = T)
#> Processing 2530 reads, using 2 gene models 
#> 2509 paired-end reads overlap 2 gene models 
#> Filtering reads by gene model... 
#> Processing sequences. This may take several minutes... 
#> Output contains: 1676 unique reads in 2 gene models

# Filter out transcripts longer than 500 bases
txReads <- tx_filter_max_width(txReads, 500) 
```

Now `txReads` contains a list with all paired-end RNA-seq mappings
divided by their corresponding gene models, along their sequences as
specified in the call to `tx_reads()`.

The resulting object is a GenomicRangesList (*GRangesList*) a list that
contains *GenomicRanges* the core object of their homonymous package,
but their coordinates belong now to the transcriptomic references used.

In this way we can take advantage of **GenomicRanges** functions and
operators to retrieve information from the mappings.

For example:

  - The names of the gene models in the list

<!-- end list -->

``` r
names(txReads)
#> [1] "uc003lam.1" "uc010nap.1"
```

  - The ranges of the mappings for an specific gene

<!-- end list -->

``` r
txReads$uc003lam.1@ranges # Ranges of reads in transcriptomic space
#> IRanges object with 1622 ranges and 0 metadata columns:
#>                              start       end     width
#>                          <integer> <integer> <integer>
#>   ID38046662_GCT_CCTATAT      1813      1918       106
#>   ID28233543_TCA_CCTATAT      1784      1918       135
#>   ID28233549_TCA_CCTATAT      1784      1918       135
#>   ID13461013_GGA_CCTATAT      1781      1924       144
#>   ID16878801_GTA_CCTATAT      1781      1868        88
#>                      ...       ...       ...       ...
#>   ID53343997_AGA_CCTATAT         8       139       132
#>    ID9459288_TGC_CCTATAT         8       139       132
#>   ID43960170_ACC_CCTATAT         4       139       136
#>   ID10702999_CCC_CCTATAT         3       139       137
#>   ID34387939_TTA_CCTATAT         1       139       139
```

  - Sequences from an individual mapping in the gene

<!-- end list -->

``` r
txReads$uc010nap.1$seq[6] 
#> [1] "ACAAGGATGGAAGAGGCCCTCGGGCCTGACAACACGC.............ATTGCCACCTACTTCGTGGCATCTAACCATCGTTTTT"
```

### Raw Gene Counts

A common task in RNA-seq analysis workflows is simply counting the reads
(or mappings) that fall into a gene model. This can be done using the
`tx_counts()` function.

``` r
tx_counts(txReads)
#> .
#> uc003lam.1 uc010nap.1 
#>       1622         54
```

### RNA-seq reads summary tables

Another useful representation of RNA-seq information is to summarise
reads metrics into tables spanning the whole transcript with information
per nucleotide. The `tx_covNucFreqDT()` function creates such tables
adding both the genomic and transcriptomic coordinates system and the
following metrics:

  - Coverage
  - Starts or 5’-ends counts
  - Ends or 3’-ends counts
  - Nucleotide frequencies
  - Deletion frequencies

<!-- end list -->

``` r
resTab1 <- tx_coverageDT(txReads, geneAnnot)
resTab1[[1]]
#>        chr   gencoor strand       gene txcoor cov start_5p end_3p
#>    1: chr5 134734928      - uc003lam.1      1   1        1      0
#>    2: chr5 134734927      - uc003lam.1      2   1        0      0
#>    3: chr5 134734926      - uc003lam.1      3   2        1      0
#>    4: chr5 134734925      - uc003lam.1      4   3        1      0
#>    5: chr5 134734924      - uc003lam.1      5   3        0      0
#>   ---                                                            
#> 1920: chr5 134670075      - uc003lam.1   1920   7        0      0
#> 1921: chr5 134670074      - uc003lam.1   1921   7        0      0
#> 1922: chr5 134670073      - uc003lam.1   1922   7        0      0
#> 1923: chr5 134670072      - uc003lam.1   1923   7        0      1
#> 1924: chr5 134670071      - uc003lam.1   1924   6        0      6

resTab2 <- tx_nucFreqDT(txReads, geneAnnot)
resTab2[[1]]
#>        chr   gencoor strand       gene txcoor A C G T N - .
#>    1: chr5 134734928      - uc003lam.1      1 1 0 0 0 0 0 0
#>    2: chr5 134734927      - uc003lam.1      2 0 1 0 0 0 0 0
#>    3: chr5 134734926      - uc003lam.1      3 0 0 0 2 0 0 0
#>    4: chr5 134734925      - uc003lam.1      4 0 0 3 0 0 0 0
#>    5: chr5 134734924      - uc003lam.1      5 0 0 3 0 0 0 0
#>   ---                                                      
#> 1920: chr5 134670075      - uc003lam.1   1920 0 7 0 0 0 0 0
#> 1921: chr5 134670074      - uc003lam.1   1921 0 7 0 0 0 0 0
#> 1922: chr5 134670073      - uc003lam.1   1922 0 0 7 0 0 0 0
#> 1923: chr5 134670072      - uc003lam.1   1923 0 0 0 7 0 0 0
#> 1924: chr5 134670071      - uc003lam.1   1924 0 0 0 6 0 0 0

resTab3 <- tx_covNucFreqDT(txReads, geneAnnot)
resTab3[[1]]
#>        chr   gencoor strand       gene txcoor cov start_5p end_3p A C G T N - .
#>    1: chr5 134734928      - uc003lam.1      1   1        1      0 1 0 0 0 0 0 0
#>    2: chr5 134734927      - uc003lam.1      2   1        0      0 0 1 0 0 0 0 0
#>    3: chr5 134734926      - uc003lam.1      3   2        1      0 0 0 0 2 0 0 0
#>    4: chr5 134734925      - uc003lam.1      4   3        1      0 0 0 3 0 0 0 0
#>    5: chr5 134734924      - uc003lam.1      5   3        0      0 0 0 3 0 0 0 0
#>   ---                                                                          
#> 1920: chr5 134670075      - uc003lam.1   1920   7        0      0 0 7 0 0 0 0 0
#> 1921: chr5 134670074      - uc003lam.1   1921   7        0      0 0 7 0 0 0 0 0
#> 1922: chr5 134670073      - uc003lam.1   1922   7        0      0 0 0 7 0 0 0 0
#> 1923: chr5 134670072      - uc003lam.1   1923   7        0      1 0 0 0 7 0 0 0
#> 1924: chr5 134670071      - uc003lam.1   1924   6        0      6 0 0 0 6 0 0 0

# This will show the summarized data table for the "uc010nap.1" gene
```

The resulting object is of class `data.table`, as one can perform many
tasks using this data structure, as well as writing them to disk in a
faster way than using the base R functions.

The resulting data.table enables easy and fast access to data, ready for
manipulation and analysis, for example, creating a barplot with the
coverage column:

  - Coverage barplot

<!-- end list -->

``` r
iGene <- "uc003lam.1"
barplot(resTab3[[iGene]]$cov, main = paste(iGene, "Coverage"),
        ylab = "Counts", xlab = iGene)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

  - Nucleotide frequency barplot

<!-- end list -->

``` r
iGene <- "uc010nap.1"
barplot(t(data.frame(resTab3[[iGene]][,c("A", "T", "G", "C", "N")])),
        col = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "black"), border = "gray",
        main = paste("Nucleotide Frequency"), ylab = "Counts", xlab = iGene)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

### Adding the reference sequence

When working with transcriptomic data one would like to easily extract
the relevant sequence to work with. To do this while working with a
summarized data.table, as the ones created above, we can use the
`tx_addRefSeq()` function.

To use the `tx_addRefSeq()` function we need a reference genome. The
BSgenome project has prepackaged several genomes for easy installation.
In this case we will use the packaged reference genome for human
“BSgenome.Hsapiens.UCSC.hg19”.

``` r
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19") # Uncomment if you need to install
# Loading reference genome
genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

# Adding the reference sequence to the summarized DT object
newDT <- tx_addRefSeqDT(resTab3, genome, geneAnnot)

# We can check how our nucleotide frequency information coincides with the
# reference sequence
newDT[[2]]
#>       chr   gencoor strand       gene txcoor refSeq cov start_5p end_3p A C  G
#>   1: chr9 137029686      - uc010nap.1      1      G  38       38      0 0 0 38
#>   2: chr9 137029685      - uc010nap.1      2      T  38        0      0 0 0  0
#>   3: chr9 137029684      - uc010nap.1      3      G  38        0      0 0 0 38
#>   4: chr9 137029683      - uc010nap.1      4      T  39        1      0 0 0  0
#>   5: chr9 137029682      - uc010nap.1      5      T  40        1      0 0 0  0
#>  ---                                                                          
#> 121: chr9 137029566      - uc010nap.1    121      T   5        0      0 0 0  0
#> 122: chr9 137029565      - uc010nap.1    122      T   5        0      0 0 0  0
#> 123: chr9 137029564      - uc010nap.1    123      T   5        0      1 0 0  0
#> 124: chr9 137029563      - uc010nap.1    124      T   4        0      0 0 0  0
#> 125: chr9 137029562      - uc010nap.1    125      T   4        0      4 0 0  0
#>       T N - .
#>   1:  0 0 0 0
#>   2: 38 0 0 0
#>   3:  0 0 0 0
#>   4: 39 0 0 0
#>   5: 40 0 0 0
#>  ---         
#> 121:  5 0 0 0
#> 122:  5 0 0 0
#> 123:  5 0 0 0
#> 124:  4 0 0 0
#> 125:  4 0 0 0
```

### Writing individual DTs to files

Additionally, storing the tables in a file for later use can be done
using the **data.table** package, which allows for fast writing of data
tables using the function `fwrite()`

``` r
# Writes datatable to file
mergedTable <- do.call(resTab3, what = rbind)
data.table::fwrite(mergedTable, "tableName.txt", sep = "\t")
```

<!-- Looking for the gene in the UCSC genome browser we can check that indeed our  -->

<!-- reconstructed sequence is that of the genome. -->

<!-- ![](https://user-images.githubusercontent.com/9357097/74606095-3c1d5480-50d6-11ea-8d01-d43e1c4ac997.png) -->

<!-- ## Graphical representations  -->

<!-- A simple graphical representation can be done to represent coverage using the nucleotide  -->

<!-- frequency tables -->

## Soon to come features:

  - **Complete function documentation & manual**
  - **How to use guide and practical cases**
  - **Graphical representations functions**
      - Coverage
      - Nucleotide frequency
      - Gene model ideogram
  - **Whole Transcript reconstruction**
  - **Meta-gene analysis tools**
      - Meta-gene plots

-----

## Current limitations:

  - Strand specific RNA-seq library preparation: State of the art
    RNA-seq protocols provide strand awareness from the original RNA.
    Currently txtools is designed for such libraries in mind, but future
    improvements will also enable processing of RNA-seq libraries which
    are not strad-aware.

  - Insertions: txtools is not able to deal with insertions. This is
    mainly because insertions are not part of the original
    trasncriptomic reference space as they would alter the length of the
    gene model. This could be fixed in future versions but is not a
    priority.

## Session Info

``` r
utils::sessionInfo()
#> R version 3.6.3 (2020-02-29)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 10240)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=English_United States.1252 
#> [2] LC_CTYPE=English_United States.1252   
#> [3] LC_MONETARY=English_United States.1252
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=English_United States.1252    
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] txtools_0.0.0.1
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.4                        pillar_1.4.3                     
#>  [3] compiler_3.6.3                    GenomeInfoDb_1.22.1              
#>  [5] XVector_0.26.0                    bitops_1.0-6                     
#>  [7] tools_3.6.3                       zlibbioc_1.32.0                  
#>  [9] digest_0.6.25                     BSgenome_1.54.0                  
#> [11] lifecycle_0.2.0                   tibble_3.0.0                     
#> [13] evaluate_0.14                     lattice_0.20-38                  
#> [15] pkgconfig_2.0.3                   rlang_0.4.5                      
#> [17] Matrix_1.2-18                     cli_2.0.2                        
#> [19] DelayedArray_0.12.2               yaml_2.2.1                       
#> [21] parallel_3.6.3                    xfun_0.12                        
#> [23] GenomeInfoDbData_1.2.2            rtracklayer_1.46.0               
#> [25] stringr_1.4.0                     dplyr_0.8.5                      
#> [27] knitr_1.28                        vctrs_0.2.4                      
#> [29] Biostrings_2.54.0                 plyranges_1.6.10                 
#> [31] S4Vectors_0.24.3                  IRanges_2.20.2                   
#> [33] tidyselect_1.0.0                  stats4_3.6.3                     
#> [35] grid_3.6.3                        data.table_1.12.8                
#> [37] glue_1.3.2                        Biobase_2.46.0                   
#> [39] R6_2.4.1                          BSgenome.Hsapiens.UCSC.hg19_1.4.0
#> [41] fansi_0.4.1                       XML_3.99-0.3                     
#> [43] BiocParallel_1.20.1               rmarkdown_2.1                    
#> [45] purrr_0.3.3                       magrittr_1.5                     
#> [47] ellipsis_0.3.0                    Rsamtools_2.2.3                  
#> [49] htmltools_0.4.0                   matrixStats_0.56.0               
#> [51] BiocGenerics_0.32.0               GenomicRanges_1.38.0             
#> [53] GenomicAlignments_1.22.1          assertthat_0.2.1                 
#> [55] SummarizedExperiment_1.16.1       stringi_1.4.6                    
#> [57] RCurl_1.98-1.1                    crayon_1.3.4
```
