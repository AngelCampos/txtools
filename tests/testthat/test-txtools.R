# Packages
library(testthat)
library(txtools)
library(GenomicAlignments)

# Loading demo data
bamFile <- system.file("extdata", "example_hg19.bam", package = "txtools")
bedFile <- system.file("extdata", "twoUCSCgenes_hg19.bed", package = "txtools")
reads <- txtools::tx_load_bam(bamFile, loadSeq = T, verbose = F,
                              yieldSize = 1000, pairedEnd = T)
geneAnnot <- txtools::tx_load_bed(bedFile) # plyranges read_bed function
txReads <- txtools::tx_reads(reads, geneAnnot, withSeq = T, verbose = F) %>%
    suppressWarnings()
unAssigned_demo <- tx_getUnassignedAlignments()
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
                     verbose = FALSE) %>% suppressWarnings()

DT <- tx_makeDT_covNucFreq(reads_PE, geneAnnot = dm3_geneAnnot, genome = dm3_genome)
DT <- tx_add_diffNucToRefRatio(DT, addDiffandTotalCols = TRUE)

# Tests
testthat::expect_equal(as.character(class(reads_PE[[1]])), "GRanges")

# tx_plot_staEndCov(DT, gene = "NM_001258475", txRange = window_around(200, 10))

# Test sk1 yeast data

assigned_aligns <- tx_reads(bam_sk1, gA_sk1, minReads = 1, withSeq = TRUE, verbose = FALSE) %>%
    suppressWarnings()

# Unassigned alignments
unAssigned <- tx_getUnassignedAlignments()
testthat::expect_equal(length(unAssigned),  103L)
tx_flushUnassigned()
testthat::expect_equal(tx_getUnassignedAlignments(), NULL)

# # Dissect unassigned
# # Overlapping
uA_o <- intersect(GenomicAlignments::findOverlaps(unAssigned@first, gA_sk1)@from,
                  GenomicAlignments::findOverlaps(GenomicAlignments::invertStrand(unAssigned@last),
                                gA_sk1)@from)
testthat::expect_equal(diff(uA_o) %>% magrittr::equals(1) %>% all(), TRUE)
testthat::expect_equal(length(uA_o), 103L) #103 unassigned and all overlapping

# All have out of exons
byGenes <- hlpr_splitReadsByGenes(unAssigned, gA_sk1, "any", 1)
anyRemain <- lapply(names(byGenes), function(iGene){
    r1_sta <- GenomicAlignments::start(unAssigned[byGenes[[iGene]]]@first)
    r1_end <- GenomicAlignments::end(unAssigned[byGenes[[iGene]]]@first)
    r2_sta <- GenomicAlignments::start(unAssigned[byGenes[[iGene]]]@last)
    r2_end <- GenomicAlignments::end(unAssigned[byGenes[[iGene]]]@last)
    exonCoors <- exonBlockGen(iGene, gA_sk1)
    which(r1_sta %in% exonCoors & r1_end %in% exonCoors & r2_sta %in% exonCoors &
        r2_end %in% exonCoors)
}) %>% unlist()
testthat::expect_equivalent(anyRemain, integer())

# Would be useful to know which part of the read falls outside of the exon structure

# # Flush unassigned
#
# tmpDT <- tx_makeDT_covNucFreq(assigned_aligns, gA_sk1, genome = genome)
# tx_plot_nucFreq(tmpDT, "YDR424C", removeInsert = T, bar_border = F)
# tx_plot_staEndCov(tmpDT, "YDR424C")
# tx_plot_nucFreq(tmpDT, "YER074W-A", removeInsert = T, bar_border = F)
# tx_plot_staEndCov(tmpDT, "YER074W-A")
#
#
