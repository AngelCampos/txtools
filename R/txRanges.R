
#' Generate exon coordinates block
#'
#' @param iGene Gene
#' @param geneAnnot Gene annotation
#'
#' @return
#' @export
#'
#' @examples
exonBlockGen <- function(iGene, geneAnnot){
  iStart <- geneAnnot[iGene,]$start
  iEnd <- geneAnnot[iGene,]$end
  iStrand <- geneAnnot[iGene,]$strand
  tmpA <- geneAnnot[iGene,]$blockStarts %>% strsplit(split = ",") %>%
    unlist() %>% as.numeric()
  tmpB <- geneAnnot[iGene,]$blockSizes %>% strsplit(split = ",") %>%
    unlist() %>% as.numeric()
  iBlocks <- sapply(1:length(tmpA), function(i){
    c(tmpA[i],(tmpA[i] + tmpB[i] -1))
  }) %>% t
  iBlocks <- iStart + iBlocks
  if(iEnd %in% as.vector(iBlocks)){
    if(iStrand == "+"){
      sapply(1:nrow(iBlocks), function(k){
        iBlocks[k,1]:iBlocks[k,2]
      }) %>% unlist %>% as.numeric
    }else if(iStrand == "-"){
      sapply(1:nrow(iBlocks), function(k){
        iBlocks[k,1]:iBlocks[k,2]
      }) %>% unlist %>% rev %>% as.numeric
    }
  }else{stop(paste("Malformed exon structure at gene", iGene))}
}

#
#' Transcriptomic range to genomic coordinates
#'
#' From transcriptomic range to genomic coordinates in BED format
#'
#' @param Tgene The gene
#' @param TrangeStart txRange start
#' @param TrangeEnd txRange end
#' @param geneAnnot The gene annotation in BED style
#'
#' @return vector
#' @export
#'
#' @examples
#' # Need to add examples
txRangeToBED <- function(Tgene, TrangeStart, TrangeEnd, geneAnnot){
    genCoor <- exonBlockGen(Tgene, geneAnnot)[TrangeStart:TrangeEnd]
    if(geneAnnot[Tgene,]$strand == "-"){
        genCoor <- rev(genCoor)
    }
    gStart <- genCoor[1]
    gEnd <- genCoor[length(genCoor)]
    blockCount <- sum(diff(genCoor) !=1) + 1
    if(sum(diff(genCoor) < 1) > 0){stop("Next exon position cannot be negative")}
    blockSizes <- paste(diff(which(c(2, diff(genCoor), 2) !=1)), collapse = ",")
    blockStarts <- paste(c(0, genCoor[which(diff(genCoor) != 1) +1] - genCoor[1]), collapse = ",")
    c(geneAnnot[Tgene,]$chr, gStart, gEnd, Tgene, 0,
      geneAnnot[Tgene,]$strand, gStart, gEnd, "0", blockCount,
      blockSizes, blockStarts)
}

#' Catenate factors example function
#'
#' This is an example function
#'
#' @param a factors A
#' @param b factors B
#'
#' @return
#' @export
#'
#' @examples
#' fbind(iris$Species[c(1, 51, 101)], PlantGrowth$group[c(1, 11, 21)])
fbind <- function(a, b) {
  factor(c(as.character(a), as.character(b)))
}
