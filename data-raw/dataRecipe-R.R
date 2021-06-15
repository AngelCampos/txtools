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
bam_sk1 <- tx_load_bam("~/WORKSPACE_wexac/m6A/SchwartzCell_data/bam/SchwartzCellWTRep1_Input_AlignedReads.bam", pairedEnd = T, loadSeq = T)
gA_sk1 <- tx_load_bed("~/BIGDATA/gene_annotations/yeast/sk1_BROAD_schragi/geneAnnot_sk1_table.bed")
# bam_sk1_mazterSeq <- tx_load_bam("~/WORKSPACE_wexac/m6A/Lib298_yeast_IME4_mazF/bam/298MazFExp3_966rep1_IP_Aligned.out.sorted.bam", pairedEnd = T, loadSeq = T)
gA_sk1 <- gA_sk1[which(gA_sk1$name == "YDR424C" | gA_sk1$name == "YER074W-A" )]
ov1 <- GenomicAlignments::findOverlaps(bam_sk1@first, gA_sk1)
ov2 <- GenomicAlignments::findOverlaps(GenomicAlignments::invertStrand(bam_sk1@last), gA_sk1)
bam_sk1 <- bam_sk1[intersect(ov2@from, ov1@from)]
genome_sk1 <- tx_load_genome("/home/labs/schwartzlab/schwartz/data/genomes/Yeast_sk1/STAR/sk1.fa")[unique(as.character(seqnames(gA_sk1)))]
use_data(bam_sk1)
use_data(gA_sk1)
use_data(genome_sk1)
