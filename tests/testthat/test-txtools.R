# Packages
library(testthat)
library(txtools)

# Loading demo data
bamFile <- system.file("extdata", "example_hg19.bam", package = "txtools")
bedFile <- system.file("extdata", "twoUCSCgenes_hg19.bed", package = "txtools")
reads <- txtools::tx_load_bam(bamFile, loadSeq = T, verbose = F,
                              yieldSize = 1000, pairedEnd = T)
geneAnnot <- txtools::tx_load_bed(bedFile) # plyranges read_bed function
txReads <- txtools::tx_reads(reads, geneAnnot, withSeq = T, verbose = F)
DTL <- txtools::tx_makeDT_covNucFreq(txReads, geneAnnot)

# Tests ########################################################################
# Length of txReads
testthat::expect_identical(length(txReads), 2L)
# Class of txReads contents
testthat::expect_identical(as.character(class(txReads[[1]])), "GRanges")
# Class data.table
testthat::expect_identical(class(DTL), c("data.table", "data.frame"))
#Check for txcoor integer continuity
testthat::expect_equivalent(DTL[DTL$gene == "uc003lam.1",]$txcoor, 1:1924)
# Test for strand of gene
testthat::expect_identical(unique(DTL[DTL$gene == "uc003lam.1",]$strand), as.factor("-"))
testthat::expect_identical(unique(DTL[DTL$gene == "uc010nap.1",]$strand), as.factor("-"))
