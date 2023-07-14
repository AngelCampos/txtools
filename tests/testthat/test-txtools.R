# Packages
library(testthat)
library(txtools)
library(GenomicAlignments)
library(magrittr)

# Loading demo data
NCORES <- 1
bamFile <- system.file("extdata", "example_hg19.bam", package = "txtools")
bedFile <- system.file("extdata", "twoUCSCgenes_hg19.bed", package = "txtools")
geneAnnot <- txtools::tx_load_bed(bedFile)
reads <- txtools::tx_load_bam(bamFile, loadSeq = T, verbose = F,
                              yieldSize = 1000, pairedEnd = T)
txReads <- txtools::tx_reads(reads, geneAnnot, withSeq = T, verbose = F) %>%
    suppressWarnings()
unAssigned_demo <- tx_get_unassignedAlignments()
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
                     nCores = NCORES,
                     minReads = 1,
                     verbose = FALSE) %>% suppressWarnings()

DT <- tx_makeDT_covNucFreq(reads_PE, geneAnnot = dm3_geneAnnot, genome = dm3_genome)
DT <- tx_add_misincRate(DT, addCounts = TRUE, minNucReads = 10)
DT <- tx_add_geneRegion(DT, dm3_geneAnnot, NCORES)
# Tests
testthat::expect_equal(as.character(class(reads_PE[[1]])), "GRanges")

# tx_plot_staEndCov(DT, gene = "NM_001258475", txRange = window_around(200, 10))

# Test sk1 yeast data
assigned_aligns <- tx_reads(bam_sk1, gA_sk1, minReads = 1, withSeq = TRUE, verbose = FALSE) %>%
    suppressWarnings()

# Retrieving and flushing unassigned alignments
unAssigned <- tx_get_unassignedAlignments()
testthat::expect_equal(length(unAssigned),  104L)
tx_flushUnassigned()
testthat::expect_equal(tx_get_unassignedAlignments(), NULL)

# # Dissect unassigned
# # Overlapping
uA_o <- intersect(GenomicAlignments::findOverlaps(unAssigned@first, gA_sk1)@from,
                  GenomicAlignments::findOverlaps(GenomicAlignments::invertStrand(unAssigned@last),
                                gA_sk1)@from)
testthat::expect_equal(diff(uA_o) %>% magrittr::equals(1) %>% all(), TRUE) ## All overlap gene annotation

# Most reads have starts or end outside of exon limits
byGenes <- hlpr_splitReadsByGenes(unAssigned, gA_sk1, "any", 1)
iGene <- names(byGenes)[1]
anyRemain <- lapply(names(byGenes), function(iGene){
    r1_sta <- GenomicAlignments::start(unAssigned[byGenes[[iGene]]]@first)
    r1_end <- GenomicAlignments::end(unAssigned[byGenes[[iGene]]]@first)
    r2_sta <- GenomicAlignments::start(unAssigned[byGenes[[iGene]]]@last)
    r2_end <- GenomicAlignments::end(unAssigned[byGenes[[iGene]]]@last)
    exonCoors <- exonBlockGen(iGene, gA_sk1)
    which(r1_sta %in% exonCoors & r1_end %in% exonCoors & r2_sta %in% exonCoors &
        r2_end %in% exonCoors)
}) %>% unlist()
testthat::expect_equivalent(anyRemain, 33L)

# Gapped alignment, doesn't match exon structure
gene_i <- "YDR424C"
remainAlign <- unAssigned[byGenes[[gene_i]]][anyRemain]
which_N <- GenomicAlignments::njunc(remainAlign@first) > 0 |
    GenomicAlignments::njunc(remainAlign@last) > 0
selAligns_r1 <- remainAlign@first[which_N]
selAligns_r2 <- remainAlign@last[which_N]
selExons <- exonGRanges(gA_sk1[gA_sk1$name == gene_i])
tmpRanges_r1 <- GenomicAlignments::cigarRangesAlongReferenceSpace(
    GenomicAlignments::cigar(selAligns_r1), ops = "M", with.ops = TRUE,
    pos = GenomicRanges::start(selAligns_r1))
tmpRanges_r2 <- GenomicAlignments::cigarRangesAlongReferenceSpace(
    GenomicAlignments::cigar(selAligns_r2), ops = "M", with.ops = TRUE,
    pos = GenomicRanges::start(selAligns_r2))
tmp1 <- S4Vectors::`%in%`(GenomicRanges::end(tmpRanges_r1), GenomicRanges::end(selExons)) |
    S4Vectors::`%in%`(GenomicRanges::start(tmpRanges_r1), GenomicRanges::start(selExons))
tmp1[unlist(lapply(tmp1, "length")) == 1] <- TRUE
tmp2 <- S4Vectors::`%in%`(GenomicRanges::end(tmpRanges_r2), GenomicRanges::end(selExons)) |
    S4Vectors::`%in%`(GenomicRanges::start(tmpRanges_r2), GenomicRanges::start(selExons))
tmp2[unlist(lapply(tmp2, "length")) == 1] <- TRUE
which_N[which_N] <- !(all(tmp1) & all(tmp2))
testthat::expect_equivalent(sum(which_N), 1) # The remaining alginment was left
#                                              out due to its gaps not meeting
#                                              the gene structure

# Check that sequences from reads that traverse exons are correctly reconstructed.
reads_sk1 <- tx_reads(bam_sk1, gA_sk1, minReads = 1, verbose = F) %>% suppressWarnings()
exongGRList <- exonGRanges(gA_sk1)
splicing_check <- lapply(1:length(gA_sk1), function(i){
    iGene <- gA_sk1$name[i]
    ovl <- findOverlaps(bam_sk1[names(reads_sk1[[iGene]])]@first, exongGRList[[iGene]], type = "any")
    ovl <- split(ovl@to, ovl@from) %>% lapply(length)
    ovIndex_1 <- names(ovl)[ovl %>% unlist() %>% is_greater_than(1)] %>% as.numeric
    ovl <- findOverlaps(bam_sk1[names(reads_sk1[[iGene]])]@last, invertStrand(exongGRList[[iGene]]), type = "any")
    ovl <- split(ovl@to, ovl@from) %>% lapply(length)
    ovIndex_2 <- names(ovl)[ovl %>% unlist() %>% is_greater_than(1)] %>% as.numeric
    bam_trExons <- bam_sk1[names(reads_sk1[[iGene]])][union(ovIndex_1, ovIndex_2)]
    DT <- tx_reads(bam_trExons, gA_sk1, withSeq = T, minReads = 1, verbose = F) %>%
        suppressWarnings() %>%
        tx_makeDT_covNucFreq(geneAnnot = gA_sk1, genome = genome_sk1) %>%
        tx_add_misincCount()
    sum(DT$misincCount)
}) %>% unlist
testthat::expect_equivalent(splicing_check, c(0, 0)) # There are some mismatches, got to check them on IGV.

# # Would be useful to know which part of the read falls outside of the exon structure
# GA_index <- indexAlignmentsByGenomicRegion(bam_sk1, tx_extend_UTR(gA_sk1, 20, 20))
# GA_index %>% lapply(length) %>% unlist
# GA_index <- indexAlignmentsByGenomicRegion(bam_sk1, gA_sk1)
