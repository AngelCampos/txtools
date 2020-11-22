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
testthat::expect_equal(length(txReads), 2)
# Class of txReads contents
testthat::expect_identical(as.character(class(txReads[[1]])), "GRanges")
# Class data.table
testthat::expect_identical(class(DTL), c("data.table", "data.frame"))
#Check for txcoor integer continuity
testthat::expect_equivalent(DTL[DTL$gene == "uc003lam.1",]$txcoor, 1:1924)
# Test for strand of gene
testthat::expect_identical(unique(DTL[DTL$gene == "uc003lam.1",]$strand), as.factor("-"))
testthat::expect_identical(unique(DTL[DTL$gene == "uc010nap.1",]$strand), as.factor("-"))


# Quick example code ###########

# Getting paths to files
BED_file <- tx_dm3_geneAnnot()
FASTA_file <- pasillaBamSubset::dm3_chr4()
PE_BAM_file <- pasillaBamSubset::untreated3_chr4()

# Loading gene annotation, genome, and alignments into R.
dm3_geneAnnot <- tx_load_bed(BED_file)
dm3_genome <- tx_load_genome(FASTA_file)
dm3_PEreads <- tx_load_bam(file = PE_BAM_file, pairedEnd = TRUE, loadSeq = T, verbose = FALSE)

reads_PE <- tx_reads(reads = dm3_PEreads,
                     geneAnnot = dm3_geneAnnot,
                     withSeq = T,
                     nCores = 2,
                     minReads = 1,
                     verbose = FALSE)

DT <- tx_makeDT_covNucFreq(reads_PE, geneAnnot = dm3_geneAnnot, genome = dm3_genome)
DT <- tx_add_diffNucToRefRatio(DT, addDiffandTotalCols = TRUE)

# Tests
testthat::expect_equal(as.character(class(reads_PE[[1]])), "GRanges")

# tx_plot_staEndCov(DT, gene = "NM_001258475", txRange = window_around(200, 10))
