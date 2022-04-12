## code to prepare `dataRecipe.R` dataset goes here
IUPAC_CODE_MAP_extended <- c(Biostrings::IUPAC_CODE_MAP,
                             structure("-" ,names =  "-"),
                             structure(paste0(Biostrings::IUPAC_CODE_MAP, "-"),
                                       names = names(Biostrings::IUPAC_CODE_MAP)))
use_data(IUPAC_CODE_MAP_extended)
IUPAC_code_2nucs <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "N", "-", ".")
use_data(IUPAC_code_2nucs)
IUPAC_code_simpl <- c("A", "C", "G", "T", "N", "-", ".")
use_data(IUPAC_code_simpl)

# Sk1 data example
library(pasillaBamSubset)
library(rtracklayer)
library(GenomicRanges)


## Not used.
## TODO: CHECK IF DISPOSABLE
# Bam file
# bam_sk1 <- tx_load_bam("~/WORKSPACE_wexac/m6A/SchwartzCell_data/bam/SchwartzCellWTRep1_Input_AlignedReads.bam", pairedEnd = T, loadSeq = T)
# bam_sk1_mazterSeq <- tx_load_bam("~/WORKSPACE_wexac/m6A/Lib298_yeast_IME4_mazF/bam/298MazFExp3_966rep1_IP_Aligned.out.sorted.bam", pairedEnd = T, loadSeq = T)
# bam_sk1 <- bam_sk1[intersect(ov2@from, ov1@from)]
# use_data(bam_sk1)
# # Gene annotation (BED file)
# gA_sk1 <- tx_load_bed("~/BIGDATA/gene_annotations/yeast/sk1_BROAD_schragi/geneAnnot_sk1_table.bed")
# ov1 <- GenomicAlignments::findOverlaps(bam_sk1@first, gA_sk1)
# ov2 <- GenomicAlignments::findOverlaps(GenomicAlignments::invertStrand(bam_sk1@last), gA_sk1)
# gA_sk1 <- gA_sk1[which(gA_sk1$name == "YDR424C" | gA_sk1$name == "YER074W-A" )]
#
# genome_sk1 <- tx_load_genome("/home/labs/schwartzlab/schwartz/data/genomes/Yeast_sk1/STAR/sk1.fa")
# genome_sk1[width(genome_sk1[-17:-18]) %>% which.min]
# genome_sk1 <- genome_sk1[unique(as.character(seqnames(gA_sk1)))]
# use_data(gA_sk1)
# use_data(genome_sk1)

# seqnames(gA_sk1)
# gA_sk1[]

FASTA_path <- dm3_chr4() # D. melanogaster chr4 built-in pasillaBamSubset
genome <- tx_load_genome(fastaFile = FASTA_path)
BAM_path <- untreated3_chr4() # Paired-end BAM file built-in pasillaBamSubset
genAligns <- tx_load_bam(file = BAM_path, pairedEnd = TRUE, loadSeq = TRUE)
fullGeneAnnot <- tx_load_bed(system.file("extdata", "dm3_chr4_RefSeqGenes_UCSC.bed", package = "txtools"))

tmpSplit <- hlpr_splitReadsByGenes(reads = genAligns, bedR = fullGeneAnnot, overlapType = "within", minReads = 200)
selGA <- fullGeneAnnot[fullGeneAnnot$name %in% names(tmpSplit)]

tmpSplit <- findOverlaps(selGA, selGA)
tmpSplit2 <- split(tmpSplit@from, x = tmpSplit@to)
selGA <- selGA[tmpSplit2[!duplicated(tmpSplit2)] %>% names() %>% as.numeric()]

set.seed(2022)
selGA <- sample(selGA, 10)
rtracklayer::export.bed(object = selGA, con = "./inst/extdata/toyGeneAnnot_Dmelan_chr4.bed")

DT_covNucFreq <- tx_makeDT_covNucFreq(genAligns, geneAnnot = geneAnnot, genome = genome)
tmpDT <- DT_covNucFreq %>% tx_add_pos() %>% tx_add_motifPresence(motif = "RRACH", nucPositions = 3)
tmpDT <- tmpDT[tmpDT$RRACH_motif_3]
data.frame(chr = tmpDT$chr,
           start = tmpDT$gencoor - 1,
           end = tmpDT$gencoor,
           name = tmpDT$pos,
           score = 0,
           strand = tmpDT$strand) %>% data.table::data.table() %>%
    data.table::fwrite(file = "annotDF.txt", sep = "\t", col.names = F)
annotSites_RRACH <- tx_load_bed("annotDF.txt")
invisible(file.remove("annotDF.txt"))
use_data(annotSites_RRACH)
rm(tmpDT)
