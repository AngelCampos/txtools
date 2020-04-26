# Tests
library(testthat)

# Loading demo data
library(txtools)
bamFile <- system.file("extdata", "example_hg19.bam", package = "txtools")
bedFile <- system.file("extdata", "twoUCSCgenes_hg19.bed", package = "txtools")
reads <- txtools::tx_load_bam(bamFile, loadSeq = T, verbose = T, yieldSize = 10000)
geneAnnot <- txtools::tx_load_bed(bedFile) # plyranges read_bed function
txReads <- txtools::tx_reads(reads, geneAnnot, withSeq = T, verbose = T)
# Length
testthat::expect_identical(length(txReads), 2L)
# Resulting class
testthat::expect_identical(as.character(class(txReads)), "CompressedGRangesList")
