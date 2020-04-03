tx_genCoorTab(txReads, geneAnnot)


genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
genome
tx_addRefSeqDT <- function(DT, fastaGenome){
    if(class(DT) == "list"){
        if(sum(class(DT[[1]]) == "data.table") == 0){
            stop("x must be of data.table class")
        }

    }
}
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
)
