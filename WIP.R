tx_genCoorTab(txReads, geneAnnot)

genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
fastaGenome <- genome
DT <- resTab3





tx_addRefSeqDT(resTab3, genome, geneAnnot)


tmpA <- geneAnnot[gene,]$blockStarts %>% strsplit(split = ",") %>% unlist() %>% as.numeric()
tmpB <- geneAnnot[gene,]$blockSizes  %>% strsplit(split = ",") %>% unlist() %>% as.numeric()
iBlocks <- (sapply(1:length(tmpA), function(i){
    c(tmpA[i],(tmpA[i] + tmpB[i] -1))
}) %>% t) + geneAnnot[gene,]$start

getGenSeqInGenome <- function(gene, geneAnnot, genome){
    tmpA <- geneAnnot[gene,]$blockStarts %>% strsplit(split = ",") %>% unlist() %>% as.numeric()
    tmpB <- geneAnnot[gene,]$blockSizes  %>% strsplit(split = ",") %>% unlist() %>% as.numeric()
    iBlocks <- (sapply(1:length(tmpA), function(i){
        c(tmpA[i],(tmpA[i] + tmpB[i] -1))
    }) %>% t) + geneAnnot[gene,]$start
    if(geneAnnot[gene,]$end %in% as.vector(iBlocks)){
        if(geneAnnot[gene,]$strand == "+"){
            sapply(1:nrow(iBlocks), function(k){
                str_sub(genome[[geneAnnot[gene,]$chr]], iBlocks[k, 1], iBlocks[k,2])
            }) %>% paste(collapse = "") %>% DNAString()
        }else if(geneAnnot[gene,]$strand == "-"){
            sapply(1:nrow(iBlocks), function(k){
                str_sub(genome[[geneAnnot[gene,]$chr]], iBlocks[k, 1], iBlocks[k,2])
            }) %>% paste(collapse = "") %>% DNAString() %>% reverseComplement()
        }
    }else{stop(paste("Malformed exon structure at gene", gene))}
}


i <- 1
resTab <- tx_covNucFreqDT(txReads, geneAnnot)
DT <- resTab$uc003lam.1
sapply(DT, class)
DT[, chr:= as.factor(chr)]
DT[, gencoor:= as.integer(gencoor)]
DT[, strand:= as.factor(strand)]
DT[, gene:= as.factor(gene)]
DT[, txcoor:= as.integer(txcoor)]

oSize(DT)
oSize(resTab$uc003lam.1)
str(resTab$uc003lam.1)
str(DT)

names(DT$uc003lam.1)
for(i in names(DT)){
    DT[, strand:= as.factor(strand)]
}
DT$uc010nap.1 %>% str


do.call(resTab, what = rbind) %>% View
