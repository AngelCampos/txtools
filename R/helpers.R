# Representing gene models as GRanges from imported BEDFILE
exonGRanges <- function(geneAnnot_GR){
    iChr <- GenomicAlignments::seqnames(geneAnnot_GR) %>% as.character()
    iStart <- GenomicRanges::start(geneAnnot_GR)
    iEnd <- GenomicRanges::end(geneAnnot_GR)
    iStrand <- GenomicRanges::strand(geneAnnot_GR)
    tmpA <- GenomicRanges::start(geneAnnot_GR$blocks) %>%
        magrittr::subtract(1) %>% magrittr::add(iStart)
    tmpB <- GenomicAlignments::width(geneAnnot_GR$blocks) %>%
        magrittr::subtract(1) %>% magrittr::add(tmpA)
    listLen <- unlist(lapply(tmpA, length))
    group <- rep(1:length(geneAnnot_GR), times = listLen)
    data.frame(rep(iChr, times = listLen),
               unlist(tmpA),
               unlist(tmpB),
               rep(iStrand, times = listLen)) %>%
        magrittr::set_colnames(c("seqnames", "start", "end", "strand")) %>%
        plyranges::as_granges() %>%
        split(group) %>%
        GenomicRanges::GRangesList() %>% magrittr::set_names(geneAnnot_GR$name)
}

# Generate exon coordinates block
# Note: It actually takes a bit less time than exonGRanges
exonBlockGen <- function(iGene, geneAnnot_GR){
    iStart <- GenomicRanges::start(geneAnnot_GR[geneAnnot_GR$name == iGene])
    iEnd <- GenomicRanges::end(geneAnnot_GR[geneAnnot_GR$name == iGene])
    iStrand <- GenomicRanges::strand(geneAnnot_GR[geneAnnot_GR$name == iGene]) %>% as.character()
    tmpA <- unlist(GenomicRanges::start(geneAnnot_GR[geneAnnot_GR$name == iGene]$blocks)) -1
    tmpB <- unlist(GenomicAlignments::width(geneAnnot_GR[geneAnnot_GR$name == iGene]$blocks))
    iBlocks <- sapply(1:length(tmpA), function(i){
        c(tmpA[i], (tmpA[i] + tmpB[i] -1))
    }) %>% t %>% magrittr::add(iStart)
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

# Reads overlapping gene models
hlpr_splitReadsByGenes <- function(reads, bedR, overlapType, minReads){
    #IgnoreStrand feature missing
    allOver_1 <- GenomicRanges::findOverlaps(reads@first, bedR, type = overlapType)
    allOver_2 <- GenomicRanges::findOverlaps(GenomicAlignments::invertStrand(reads@last),
                                             bedR,type = overlapType)
    split_1 <- split(allOver_1@from, allOver_1@to)
    split_2 <- split(allOver_2@from, allOver_2@to)
    names(split_1) <- bedR[names(split_1) %>% as.numeric()]$name
    names(split_2) <- bedR[names(split_2) %>% as.numeric()]$name
    inBoth <- intersect(names(split_1), names(split_2))
    split_1 <- split_1[inBoth]
    split_2 <- split_2[inBoth]
    split_i <- vIntersect(split_1, split_2)
    names(split_i) <- names(split_1)
    split_i <- split_i[unlist(lapply(split_i, length)) %>%
                           magrittr::is_weakly_greater_than(minReads) %>% which]
    return(split_i)
}

# For single-end reads
hlpr_splitReadsByGenes_singleEnd <- function(reads, bedR, overlapType, minReads){
    allOver_1 <- GenomicRanges::findOverlaps(reads, bedR, type = overlapType)
    split_1 <- split(allOver_1@from, allOver_1@to)
    names(split_1) <- bedR[as.numeric(names(split_1))]$name
    split_1 <- split_1[sapply(split_1, length) %>%
                           magrittr::is_weakly_greater_than(minReads) %>%
                           which]
    return(split_1)
}

# Process reads into transcripts
hlpr_ReadsInGene <- function(reads, iGene, geneAnnot, split_i, allExons, minReads, withSeq){
    iStrand <- geneAnnot[which(geneAnnot$name == iGene)] %>%
        GenomicRanges::strand() %>% as.character()
    iExon <- exonBlockGen(iGene = iGene, geneAnnot_GR = geneAnnot)
    selReadsbyPair <- split_i[[iGene]]
    # Selecting paired reads to merge
    iReads_r1 <- reads@first[selReadsbyPair]
    iReads_r2 <- reads@last[selReadsbyPair]
    # Filtering reads strictly inside exons
    # Both ends of both reads fall into the gene model
    stEndTable <- as.matrix(data.frame(r1_S = GenomicRanges::start(iReads_r1),
                                       r1_E = GenomicRanges::end(iReads_r1),
                                       r2_S = GenomicRanges::start(iReads_r2),
                                       r2_E = GenomicRanges::end(iReads_r2)))
    pass <- (stEndTable %in% iExon) %>% matrix(ncol = 4, byrow = F) %>%
        rowSums %>% magrittr::equals(4) %>% which()
    if(length(pass) < minReads){return(GenomicRanges::GRanges())} # No reads Return empty GA
    # Reads cover consecutive exons
    tmp <- GenomicRanges::findOverlaps(iReads_r1[pass], allExons[[iGene]])
    passPos <- split(tmp@to, tmp@from) %>% lapply(diff) %>%
        lapply(function(x) all(x == 1)) %>% unlist %>% which
    tmp <- GenomicRanges::findOverlaps(GenomicAlignments::invertStrand(iReads_r2[pass]),
                                       allExons[[iGene]])
    passNeg <- split(tmp@to, tmp@from) %>% lapply(diff) %>%
        lapply(function(x) all(x == 1)) %>% unlist %>% which
    pass <- pass[as.numeric(intersect(passNeg, passPos))]
    if(length(pass) < minReads){return(GenomicRanges::GRanges())} # No reads Return empty GA
    # Boundaries of merged reads
    if(iStrand == "+"){
        tReads <- data.frame(start    = match(GenomicRanges::start(iReads_r1[pass]), iExon),
                             end      = match(GenomicRanges::end(iReads_r2[pass]), iExon),
                             strand   = GenomicRanges::strand(iReads_r1[pass]),
                             seqnames = iGene)
    }else if(iStrand == "-"){
        tReads <- data.frame(start    = match(GenomicRanges::end(iReads_r1[pass]), iExon),
                             end      = match(GenomicRanges::start(iReads_r2[pass]), iExon),
                             strand   = GenomicRanges::strand(iReads_r1[pass]),
                             seqnames = iGene)
    }
    pass <- pass[tReads$end >= (tReads$start -1)]
    tReads <- check_DFforGRanges(tReads) %>% plyranges::as_granges()
    names(tReads) <- names(reads)[selReadsbyPair][pass]
    GenomeInfoDb::seqlengths(tReads) <- length(iExon)
    # Result if no sequence is input or required
    if(withSeq == F){
        return(tReads)
    }
    # Merging sequences
    tReads$start_r1 <- GenomicRanges::start(iReads_r1[pass])
    tReads$end_r1 <- GenomicRanges::end(iReads_r1[pass])
    tReads$start_r2 <- GenomicRanges::start(iReads_r2[pass])
    tReads$end_r2 <- GenomicRanges::end(iReads_r2[pass])
    tReads$strand_r1 <- GenomicRanges::strand(iReads_r1[pass])
    tReads$cigar_r1 <- GenomicAlignments::cigar(iReads_r1[pass])
    tReads$cigar_r2 <- GenomicAlignments::cigar(iReads_r2[pass])
    tReads$seq_r1 <- S4Vectors::mcols(iReads_r1[pass])$seq
    tReads$seq_r2 <- S4Vectors::mcols(iReads_r2[pass])$seq
    # Constructing the merged read sequence
    if(iStrand == "+"){
        tReads$seq1 <- GenomicAlignments::sequenceLayer(S4Vectors::mcols(tReads)$seq_r1,
                                                        S4Vectors::mcols(tReads)$cigar_r1)
        tReads$seq2 <- GenomicAlignments::sequenceLayer(S4Vectors::mcols(tReads)$seq_r2,
                                                        S4Vectors::mcols(tReads)$cigar_r2)
    }else if(iStrand == "-"){
        tReads$seq1 <- GenomicAlignments::sequenceLayer(S4Vectors::mcols(tReads)$seq_r1,
                                                        S4Vectors::mcols(tReads)$cigar_r1) %>%
            Biostrings::reverseComplement()
        tReads$seq2 <- GenomicAlignments::sequenceLayer(S4Vectors::mcols(tReads)$seq_r2,
                                                        S4Vectors::mcols(tReads)$cigar_r2) %>%
            Biostrings::reverseComplement()
    }
    # Trim overflowing reads
    i <- which(Biostrings::nchar(tReads$seq1) > GenomicAlignments::width(tReads))
    tReads[i]$seq1 <- stringr::str_sub(tReads[i]$seq1, start = 1, GenomicAlignments::width(tReads[i]))
    i <- which(Biostrings::nchar(tReads$seq2) > GenomicAlignments::width(tReads))
    tReads[i]$seq2 <- stringr::str_sub(tReads[i]$seq2, start = - GenomicAlignments::width(tReads[i]))
    # Calculate overlap
    tReads$diffSeq <- Biostrings::nchar(tReads$seq1) + Biostrings::nchar(tReads$seq2) - GenomicAlignments::width(tReads)
    tReads$oL <- tReads$diffSeq %>% magrittr::is_greater_than(0)
    tReads$mergedSeq <- "." # Place holder
    # No overlap reads
    i <- which(!(tReads$oL))
    if(length(i) > 0){
        gapLen <- GenomicAlignments::width(tReads[i]) -
            Biostrings::nchar(tReads[i]$seq1) - Biostrings::nchar(tReads[i]$seq2)
        insSeq <- stringr::str_dup(string = ".", gapLen)
        tReads[i]$mergedSeq <- paste(tReads[i]$seq1, insSeq,
                                     tReads[i]$seq2, sep = "") %>%
            Biostrings::DNAStringSet()
    }
    # Overlapped reads
    i <- which(tReads$oL)
    if(length(i) > 0){
        tmp1 <- stringr::str_sub(tReads[i]$seq1, start = -tReads[i]$diffSeq)
        tmp2 <- stringr::str_sub(tReads[i]$seq2, start = 1, end = tReads[i]$diffSeq)
        ovSeq <- rep(".", length(i))
        for(j in 1:length(i)){
            if(tmp1[j] == tmp2[j]){
                ovSeq[j] <- tmp1[j]
            }else{
                ovSeq[j] <- Biostrings::DNAStringSet(c(tmp1[j], tmp2[j])) %>%
                    Biostrings::consensusMatrix() %>%
                    Biostrings::consensusString(ambiguityMap = txtools::IUPAC_CODE_MAP_extended,
                                                threshold = 0.16)
            }
        }
        tReads[i]$mergedSeq <- paste0(stringr::str_sub(tReads[i]$seq1,
                                                       start = 1,
                                                       end = Biostrings::nchar(tReads[i]$seq1) -
                                                           tReads[i]$diffSeq),
                                      ovSeq, stringr::str_sub(tReads[i]$seq2,
                                                              start = tReads[i]$diffSeq + 1,
                                                              end = Biostrings::nchar(tReads[i]$seq2)))
    }
    # Final reads
    fReads <- tReads
    S4Vectors::mcols(fReads) <- NULL
    fReads$seq <- tReads$mergedSeq
    return(fReads)
}

# Vectorized version of helper
v_hlpr_ReadsInGene <- Vectorize(hlpr_ReadsInGene, "iGene")

# reads in gene single-end
hlpr_ReadsInGene_SingleEnd <- function(reads, iGene, geneAnnot, split_i,
                                       allExons, minReads, withSeq){
    iStrand <- geneAnnot[which(geneAnnot$name == iGene)] %>%
        GenomicRanges::strand() %>% as.character()
    iExon <- exonBlockGen(iGene = iGene, geneAnnot_GR = geneAnnot)
    selReadsbyPair <- split_i[[iGene]]
    # Selecting paired reads to merge
    iReads_r1 <- reads[selReadsbyPair]
    # Filtering reads extrictly inside exons
    # Both ends must fall into the gene model
    stEndTable <- as.matrix(data.frame(r1_S = GenomicRanges::start(iReads_r1),
                                       r1_E = GenomicRanges::end(iReads_r1)))
    pass <- (stEndTable %in% iExon) %>% matrix(ncol = 2, byrow = F) %>%
        rowSums() %>% magrittr::equals(2) %>% which()
    if(length(pass) < minReads){return(GenomicRanges::GRanges())} # No reads Return empty GA
    # Boundaries of merged reads
    if(iStrand == "+"){
        tReads <- data.frame(start  = match(GenomicRanges::start(iReads_r1[pass]), iExon),
                             end      = match(GenomicRanges::end(iReads_r1[pass]), iExon),
                             strand   = GenomicRanges::strand(iReads_r1[pass]),
                             seqnames = iGene)
    }else if(iStrand == "-"){
        tReads <- data.frame(start    = match(GenomicRanges::end(iReads_r1[pass]), iExon),
                             end      = match(GenomicRanges::start(iReads_r1[pass]), iExon),
                             strand   = GenomicRanges::strand(iReads_r1[pass]),
                             seqnames = iGene)
    }
    pass <- pass[tReads$end >= (tReads$start -1)]
    tReads <- check_DFforGRanges(tReads) %>% plyranges::as_granges()
    names(tReads) <- names(reads)[selReadsbyPair][pass]
    GenomeInfoDb::seqlengths(tReads) <- length(iExon)
    # Removing reads that don't match in length when passed to transcriptomic space
    R1_len <- GenomicAlignments::cigarWidthAlongReferenceSpace(
        GenomicAlignments::cigar(iReads_r1[pass]), N.regions.removed = T)
    tReads_c <- tReads[which(GenomicRanges::width(tReads) == R1_len)]
    pass <- pass[which(GenomicRanges::width(tReads) == R1_len)]
    if(length(pass) < minReads){return(GenomicRanges::GRanges())} # No reads Return empty GA
    tReads <- tReads_c
    # Result if no sequence is input or required
    if(withSeq == F){
        return(tReads)
    }
    # Merging sequences
    tReads$start_r1 <- GenomicRanges::start(iReads_r1[pass])
    tReads$end_r1 <- GenomicRanges::end(iReads_r1[pass])
    tReads$strand_r1 <- GenomicRanges::strand(iReads_r1[pass])
    tReads$cigar_r1 <- GenomicAlignments::cigar(iReads_r1[pass])
    tReads$seq_r1 <- S4Vectors::mcols(iReads_r1[pass])$seq
    # Constructing the merged read sequence
    if(iStrand == "+"){
        tReads$seq1 <- GenomicAlignments::sequenceLayer(S4Vectors::mcols(tReads)$seq_r1,
                                                        S4Vectors::mcols(tReads)$cigar_r1,
                                                        to = "reference-N-regions-removed")
    }else if(iStrand == "-"){
        tReads$seq1 <- Biostrings::reverseComplement(
            GenomicAlignments::sequenceLayer(S4Vectors::mcols(tReads)$seq_r1,
                                             S4Vectors::mcols(tReads)$cigar_r1,
                                             to = "reference-N-regions-removed"))

    }
    # Trim overflowing reads
    i <- which(Biostrings::nchar(tReads$seq1) > GenomicAlignments::width(tReads))
    if(length(i) > 0){
        tReads[i]$seq1 <- stringr::str_sub(tReads[i]$seq1, start = 1,
                                           GenomicAlignments::width(tReads[i]))
    }
    # Calculate overlap
    tReads$diffSeq <- Biostrings::nchar(tReads$seq1) - GenomicAlignments::width(tReads)
    tReads$oL <- magrittr::is_greater_than(tReads$diffSeq, 0)
    tReads$mergedSeq <- "." # Place holder
    # No overlap reads
    i <- which(!(tReads$oL))
    if(length(i) > 0){
        gapLen <- GenomicAlignments::width(tReads[i]) - Biostrings::nchar(tReads[i]$seq1)
        insSeq <- stringr::str_dup(string = ".", gapLen)
        tReads[i]$mergedSeq <- paste(tReads[i]$seq1, insSeq, sep = "") %>%
            Biostrings::DNAStringSet()
    }
    # Final reads
    fReads <- tReads
    S4Vectors::mcols(fReads) <- NULL
    fReads$seq <- tReads$mergedSeq
    return(fReads)
}

# Calculate coverage table: coverage, 5prime-starts, and 3prime-ends
hlp_coverageTab <- function(x){
    if(class(x) != "CompressedGRangesList"){
        stop("x must be of class CompressedGRangesList")
    }
    if(names(x) %>% duplicated() %>% sum() %>% magrittr::is_greater_than(0)){
        stop("List contains duplicated gene models")
    }
    cov <- hlp_coverage(x) %>% lapply(as.vector)
    sta <- GenomicRanges::start(x) %>% lapply(table) %>% magrittr::set_names(names(x))
    end <- GenomicRanges::end(x) %>% lapply(table) %>% magrittr::set_names(names(x))
    len <- GenomeInfoDb::seqlengths(x)
    OUT <- lapply(names(x), function(i){
        tmpMat <- matrix(0, nrow = len[[i]], ncol = 3) %>% data.frame() %>%
            data.table::data.table() %>%
            magrittr::set_rownames(1:len[[i]]) %>%
            magrittr::set_colnames(c("cov", "start_5p", "end_3p"))
        tmpMat[names(sta[[i]]), "start_5p"] <- sta[[i]]
        tmpMat[names(end[[i]]), "end_3p"] <- end[[i]]
        tmpMat[, "cov"] <- cov[[i]] %>% as.numeric()
        tmpMat$cov <- as.integer(tmpMat$cov)
        tmpMat$start_5p <- as.integer(tmpMat$start_5p)
        tmpMat$end_3p <- as.integer(tmpMat$end_3p)
        return(tmpMat)
    }) %>% magrittr::set_names(names(x))
    return(OUT)
}

#' Calculate coverage table (coverage, ends, starts) for all gene models
#' without gene coordinates (Multicore)
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via tx_reads().
#' @param nCores integer. Number of cores to use to run function.
#'
#' @return data.table
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
hlp_coverageTab_mc <- function(x, nCores){
    if(class(x) != "CompressedGRangesList"){
        stop("x must be of class SimpleGRangesList")
    }
    if(names(x) %>% duplicated() %>% sum() %>% magrittr::is_greater_than(0)){
        stop("List contains duplicated gene models")
    }
    cov <- hlp_coverage(x) %>% parallel::mclapply(as.vector, mc.cores = nCores)
    sta <- GenomicRanges::start(x) %>% lapply(table) %>% magrittr::set_names(names(x))
    end <- GenomicRanges::end(x) %>% lapply(table) %>% magrittr::set_names(names(x))
    len <- lapply(cov, length) %>% unlist
    OUT <- parallel::mclapply(mc.cores = nCores, 1:length(x), function(i){
        tmpMat <- matrix(0, nrow = len[i], ncol = 3) %>% data.frame() %>%
            magrittr::set_rownames(1:len[i]) %>%
            magrittr::set_colnames(c("cov", "start_5p", "end_3p"))
        tmpMat[names(sta[[i]]), "start_5p"] <- sta[[i]]
        tmpMat[names(end[[i]]), "end_3p"] <- end[[i]]
        tmpMat[, "cov"] <- cov[[i]] %>% as.numeric()
        tmpMat <- tmpMat %>% data.table::data.table()
        tmpMat$cov <- as.integer(tmpMat$cov)
        tmpMat$start_5p <- as.integer(tmpMat$start_5p)
        tmpMat$end_3p <- as.integer(tmpMat$end_3p)
        return(tmpMat)
    }) %>% magrittr::set_names(names(x))
    return(OUT)
}

#' Calculate nucleotide frequency pileup for all gene models
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via tx_reads().
#' @param simplify_IUPAC character. Method to simplify ambiguous reads:
#' 1) 'not': Ambiguous reads will be left as their IUPAC_ambiguous code, e.g.
#' for positions in which a 'G' and and 'A' where read an 'R' will note this.
#' 2) "splitHalf': Ambiguous reads will be divided in half, in odd cases having
#' fractions of reads assigned.
#' 3) "splitForceInt": Ambiguous reads will be divided in half, but forced to be
#' integers, unassigned fraction of reads will be summed and assigned as 'N'.
#'
#' @return data.table
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
hlp_nucFreqTab <- function(x, simplify_IUPAC = "not"){
    if(!all(simplify_IUPAC %in% c("not", "splitHalf", "splitForceInt"))){
        stop("simplify_IUPAC argument must be either: 'not', 'splitHalf' or ",
             "'splitForceInt'.")
    }
    iGenes <- names(x)
    lapply(iGenes, function(iGene){
        y <- Biostrings::consensusMatrix(x = x[[iGene]]$seq,
                                         shift = GenomicRanges::start(x[[iGene]]) -1,
                                         width = GenomeInfoDb::seqlengths(x[[iGene]])[iGene])
        if(simplify_IUPAC == "not"){
            hlp_addMissingNucs(y) %>% t %>% data.table::data.table()
        }else if(simplify_IUPAC == "splitHalf"){
            hlp_splitNucsHalf(y) %>% t %>% data.table::data.table()
        }else if(simplify_IUPAC == "splitForceInt"){
            hlp_splitNucsForceInt(y) %>% t %>% data.table::data.table()
        }
    }) %>% magrittr::set_names(names(x))
}

#' Calculate nucleotide frequency pileup for all gene models
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via tx_reads().
#' @param simplify_IUPAC character. Method to simplify ambiguous reads:
#' 1) 'not': Ambiguous reads will be left as their IUPAC_ambiguous code, e.g.
#' for positions in which a 'G' and and 'A' where read an 'R' will note this.
#' 2) "splitHalf': Ambiguous reads will be divided in half, in odd cases having
#' fractions of reads assigned.
#' 3) "splitForceInt": Ambiguous reads will be divided in half, but forced to be
#' integers, unassigned fraction of reads will be summed and assigned as 'N'.
#' @param nCores integer. Number of cores to use to run function.
#'
#' @return data.table
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
hlp_nucFreqTab_mc <- function(x, simplify_IUPAC = "not", nCores){
    iGenes <- names(x)
    parallel::mclapply(mc.cores = nCores, X = iGenes, function(iGene){
        y <- Biostrings::consensusMatrix(x = x[[iGene]]$seq,
                                         shift = GenomicRanges::start(x[[iGene]]) -1,
                                         width = GenomeInfoDb::seqlengths(x[[iGene]])[iGene])
        if(simplify_IUPAC == "splitForceInt"){
            hlp_splitNucsForceInt(y) %>% t %>% apply(MARGIN = 2, FUN = "as.integer") %>% data.table::data.table()
        }else if(simplify_IUPAC == "splitHalf"){
            hlp_splitNucsHalf(y) %>% t %>% data.table::data.table()
        }else if(simplify_IUPAC == "not"){
            hlp_addMissingNucs(y) %>% t %>% data.table::data.table()
        }
    }) %>% magrittr::set_names(names(x))
}

#' Table with genomic and transcriptomic coordinates
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via tx_reads().
#' @param geneAnnot GenomicRanges. Gene annotation loaded via the tx_load_bed()
#'
#' @return data.table
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
hlp_genCoorTab <- function(x, geneAnnot){
    if(all(names(x) %in% geneAnnot$name)){
        lapply(names(x), function(iGene){
            tmp2 <- geneAnnot[which(geneAnnot$name == iGene)]
            tmp3 <- c(GenomicAlignments::seqnames(tmp2),
                      GenomicRanges::strand(tmp2)) %>% as.character() %>% c(iGene)
            tmpDT <- rep(tmp3, GenomeInfoDb::seqlengths(x[[iGene]])[iGene]) %>%
                matrix(ncol = 3, byrow = T) %>%
                cbind(exonBlockGen(iGene, geneAnnot)) %>%
                cbind(seq(1, GenomeInfoDb::seqlengths(x[[iGene]])[iGene]))
            tmpDT <- tmpDT[,c(1,4,2,3,5)] %>% data.table::data.table() %>%
                magrittr::set_colnames(c("chr", "gencoor", "strand", "gene", "txcoor"))
            tmpDT$chr <- as.factor(tmpDT$chr)
            tmpDT$gencoor <- as.integer(tmpDT$gencoor)
            tmpDT$strand <- as.factor(tmpDT$strand)
            tmpDT$gene <- as.factor(tmpDT$gene)
            tmpDT$txcoor <- as.integer(tmpDT$txcoor)
            return(tmpDT)
        }) %>% magrittr::set_names(names(x))
    }else{
        stop("Names of x are not contained in geneAnnot$name")
    }
}

#' Table with genomic and transcriptomic coordinates (Multi-core)
#'
#' This function is the multi-core version of hlp_genCoorTab(), and is only
#' available for use in UNIX-like operative systems.
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via tx_reads().
#' @param geneAnnot GenomicRanges. Gene annotation loaded via the tx_load_bed()
#' @param nCores integer. Number of cores to use to run function.
#'
#' @return data.table
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
hlp_genCoorTab_mc <- function(x, geneAnnot, nCores){
    check_mc_windows(nCores)
    if(all(names(x) %in% geneAnnot$name)){
        parallel::mclapply(mc.cores = nCores, X = names(x), function(iGene){
            tmp2 <- geneAnnot[which(geneAnnot$name == iGene)]
            tmp3 <- c(GenomicAlignments::seqnames(tmp2),
                      GenomicRanges::strand(tmp2)) %>% as.character() %>% c(iGene)
            tmpDT <- rep(tmp3, GenomeInfoDb::seqlengths(x[[iGene]])[iGene]) %>%
                matrix(ncol = 3, byrow = T) %>%
                cbind(exonBlockGen(iGene, geneAnnot)) %>%
                cbind(seq(1, GenomeInfoDb::seqlengths(x[[iGene]])[iGene]))
            tmpDT <- tmpDT[,c(1,4,2,3,5)] %>% data.table::data.table() %>%
                magrittr::set_colnames(c("chr", "gencoor", "strand", "gene", "txcoor"))
            tmpDT$chr <- as.factor(tmpDT$chr)
            tmpDT$gencoor <- as.integer(tmpDT$gencoor)
            tmpDT$strand <- as.factor(tmpDT$strand)
            tmpDT$gene <- as.factor(tmpDT$gene)
            tmpDT$txcoor <- as.integer(tmpDT$txcoor)
            return(tmpDT)
        }) %>% magrittr::set_names(names(x))
    }else{
        stop("Names of x are not contained in geneAnnot$name .\n")
    }
}

# Calculate coverage for all gene models
hlp_coverage <- function(x){
    suppressWarnings(unlist(x)) %>% GenomicRanges::coverage()
}

# Add missing nucleotides in nuc frequency table
hlp_addMissingNucs <- function(x){
    misNucs <- txtools::IUPAC_code_2nucs[which(!(txtools::IUPAC_code_2nucs %in% rownames(x)))]
    if(length(misNucs) > 0){
        tmp <- matrix(0, nrow = length(misNucs), ncol = ncol(x)) %>%
            magrittr::set_rownames(misNucs) %>% rbind(x)
        tmp[txtools::IUPAC_code_2nucs,]
    }else{
        x[txtools::IUPAC_code_2nucs,]
    }
}

# Split nucleotides in half
hlp_splitNucsHalf <- function(x){
    misNucs <- txtools::IUPAC_code_simpl[which(!(txtools::IUPAC_code_simpl %in% rownames(x)))]
    altNucs <- intersect(txtools::IUPAC_code_2nucs[5:10], rownames(x))
    x <- hlp_addMissingNucs(x)
    for(i in altNucs){
        resNuc <- Biostrings::IUPAC_CODE_MAP[i] %>% stringr::str_split("") %>% unlist
        x[resNuc[1],] <- x[resNuc[1],] + (x[i, ] / 2)
        x[resNuc[2],] <- x[resNuc[2],] + (x[i, ] / 2)
    }
    x[txtools::IUPAC_code_simpl,]
}

# Split nucleotides in half, forcing integers
hlp_splitNucsForceInt <- function(x){
    misNucs <- txtools::IUPAC_code_simpl[which(!(txtools::IUPAC_code_simpl %in% rownames(x)))]
    altNucs <- intersect(txtools::IUPAC_code_2nucs[5:10], rownames(x))
    x <- hlp_addMissingNucs(x)
    for(i in altNucs){
        if(all(x[i,] %% 2 == 0)){
            resNuc <- Biostrings::IUPAC_CODE_MAP[i] %>% stringr::str_split("") %>% unlist
            x[resNuc[1],] <- x[resNuc[1],] + (x[i,] / 2)
            x[resNuc[2],] <- x[resNuc[2],] + (x[i,] / 2)
        }else{
            resNuc <- Biostrings::IUPAC_CODE_MAP[i] %>% stringr::str_split("") %>% unlist
            x[resNuc[1],] <- x[resNuc[1],] + (x[i,] %>% magrittr::divide_by(2) %>% floor)
            x[resNuc[2],] <- x[resNuc[2],] + (x[i,] %>% magrittr::divide_by(2) %>% floor)
            x["N", ] <- x["N",] + (x[i,] - (x[i,] %>% magrittr::divide_by(2) %>%
                                                floor %>% magrittr::multiply_by(2)))
        }
    }
    x[txtools::IUPAC_code_simpl,]
}

# Helper bind two results tables
hlp_cbind2Tabs <- function(gencoorT, tab1){
    if(all(names(gencoorT) == names(tab1))){
        lapply(seq(1, length(gencoorT)), function(i){
            cbind(gencoorT[[i]], tab1[[i]])
        }) %>% magrittr::set_names(names(gencoorT))
    }
}

# Helper bind three results tables
hlp_cbind3Tabs <- function(gencoorT, tab1, tab2){
    if(all(names(gencoorT) == names(tab1) &
           names(tab2) == names(gencoorT))){
        lapply(seq(1, length(gencoorT)), function(i){
            cbind(gencoorT[[i]], tab1[[i]], tab2[[i]])
        }) %>% magrittr::set_names(names(gencoorT))
    }
}

# Unlist if IrangesList
if_IRangesList_Unlist <- function(x){
    if(methods::is(x, "IRangesList")){
        unlist(x)
    }else{
        x
    }
}

# Helper add reference sequence to DT
hlp_add_refSeqDT <- function (DT, genome, geneAnnot){
    iGene <- as.character(DT$gene[1])
    iChr <- as.character(DT$chr[1])
    iStr <- as.character(DT$strand[1])
    iGA <- geneAnnot[geneAnnot$name == iGene]
    iBlocks <- S4Vectors::mcols(iGA)$blocks %>% if_IRangesList_Unlist() %>%
        IRanges::shift(IRanges::start(iGA) - 1)
    tmp <- Biostrings::DNAString(
        paste(collapse = "",
              stringr::str_sub(string = genome[[iChr]],
                               start = IRanges::start(iBlocks),
                               end = IRanges::end(iBlocks))))
    if (iStr == "-") {
        tmp <- Biostrings::reverseComplement(tmp)
    }
    tmp <- stringr::str_split(as.character(tmp), "") %>% unlist()
    tibble::add_column(DT, refSeq = tmp, .after = "txcoor")
}

# Remove reads with nucleotide sequence length of 0 (bug detected when using Nanopore reads from Amit)
hlp_cleanBam_emptySeq <- function(reads, verbose){
    if(methods::is(reads, class2 = "GAlignments")){
        lSel <- BiocGenerics::width(S4Vectors::mcols(reads)$seq) != 0
        if(verbose){cat(length(reads) - sum(lSel), "reads filtered out for empty sequence field")}
        reads[lSel]
    }else if(methods::is(reads, class2 = "GAlignmentPairs")){
        lSel <- BiocGenerics::width(S4Vectors::mcols(reads@first)$seq) != 0 &
            BiocGenerics::width(S4Vectors::mcols(reads@last)$seq) != 0
        if(verbose){cat(length(reads) - sum(lSel), "reads filtered out for empty sequence field")}
        reads[lSel]
    }
}

.dump_envir_txtools <- new.env(parent=emptyenv())
.dumpEnvirTxtools <- function() .dump_envir_txtools

# Dump not assigned alignments to dump environment
# reads <- dm3_PEreads
# OUT <- reads_SE
hlp_dump_notAssigned <- function(reads, OUT){
    assign("notAssignedAlignments",
           reads[!(names(reads) %in% unique(unlist(lapply(OUT, names))))],
           envir = .dumpEnvirTxtools())
}

# Flush unassigned
tx_flushUnassigned <- function(){
    objnames <- ls(envir = .dumpEnvirTxtools())
    rm(list = objnames, envir = .dumpEnvirTxtools())
}

# Get vector elements from start to end
hlp_getVectorElements <- function(x, start, end){
    x[start:end]
}
get_VectorElements <- Vectorize(hlp_getVectorElements, c("start", "end"), SIMPLIFY = FALSE)


# tx_add_XXX() #################################################################

# Mark motif location in tx_DT as TRUE, it can be set for specific nuc positions or all
hlp_add_motifPresence <- function(DT, motif, nucPositions, midMot){
    DTseq <- paste(DT$refSeq, collapse = "")
    tmpLoc <- Biostrings::matchPattern(Biostrings::DNAString(motif),
                                       Biostrings::DNAString(DTseq), fixed = F)
    tmpLoc <- data.frame(tmpLoc@ranges)
    addMotif <- rep(FALSE, nrow(DT))
    if(is.numeric(nucPositions)){
        if(!nrow(tmpLoc) == 0){
            for(i in 1:nrow(tmpLoc)){
                addMotif[(tmpLoc[i, 1]) + nucPositions - 1] <- TRUE
            }
        }
    }else if(nucPositions == "all"){
        for(i in 1:nrow(tmpLoc)){
            addMotif[(tmpLoc[i, 1]):(tmpLoc[i, 2])] <- TRUE
        }
    }else if(nucPositions == "center"){
        for(i in 1:nrow(tmpLoc)){
            addMotif[(tmpLoc[i, 1]) + midMot - 1] <- TRUE
        }
    }
    tibble::add_column(DT, addMotif)
}

# Helper add annotation of sites in GRanges
hlp_add_siteAnnotation <- function (x, GRanges, colName){
    subGR <- GRanges[as.character(x$chr[1]) == GenomicRanges::seqnames(GRanges)]
    oNames <- names(x)
    addAnnot <- rep(FALSE, nrow(x))
    if(length(subGR) == 0){
        tibble::add_column(x, addAnnot) %>% magrittr::set_names(c(oNames, colName))
    }
    foundGenLoc <- GenomicRanges::start(subGR)[
        which(as.logical((GenomicRanges::start(subGR) %in% x$gencoor) &
                             (GenomicRanges::strand(subGR) == as.character(unique(x$strand)))))]
    if(length(foundGenLoc) == 0){
        tibble::add_column(x, addAnnot) %>% magrittr::set_names(c(oNames,
                                                                  colName))
    }else{
        addAnnot[match(foundGenLoc, x$gencoor)] <- TRUE
        tibble::add_column(x, addAnnot) %>% magrittr::set_names(c(oNames, colName))
    }
}

# Vectorized intersect
vIntersect <- Vectorize(intersect, c("x", "y"), SIMPLIFY = F)

# Object size
oSize <- function(x){
    print(utils::object.size(x), units = "auto")
}

# Stretching 5p-most blocks
stretchBlocks_5p <- function(blocks, extend, strand){
    blocks <- GenomicRanges::shift(blocks, shift = ifelse(strand == "+", extend, 0))
    part <- IRanges::PartitioningByEnd(blocks)
    collapsed <- unlist(blocks)
    ind_5p_pos <- IRanges::start(part)[as.logical(strand == "+")]
    ind_5p_neg <- IRanges::end(part)[as.logical(strand == "-")]
    IRanges::start(collapsed[ind_5p_pos]) <- IRanges::start(collapsed[ind_5p_pos]) - extend
    IRanges::end(collapsed[ind_5p_neg]) <- IRanges::end(collapsed[ind_5p_neg]) + extend
    return(IRanges::relist(collapsed, part))
}

# Stretching 3p-most blocks
stretchBlocks_3p <- function(blocks, extend, strand){
    blocks <- GenomicRanges::shift(blocks, shift = ifelse(strand == "-", extend, 0))
    part <- IRanges::PartitioningByEnd(blocks)
    collapsed <- unlist(blocks)
    ind_3p_pos <- IRanges::start(part)[as.logical(strand == "-")]
    ind_3p_neg <- IRanges::end(part)[as.logical(strand == "+")]
    IRanges::start(collapsed[ind_3p_pos]) <- IRanges::start(collapsed[ind_3p_pos]) - extend
    IRanges::end(collapsed[ind_3p_neg]) <- IRanges::end(collapsed[ind_3p_neg]) + extend
    return(IRanges::relist(collapsed, part))
}

# Generating DTs ###############################################################
# Generates transcriptomic coordinates table from a list of genes
hlpr_genCoorTabGenes <- function(genes, geneAnnot, genome = NULL, nCores = 1){
    if(all(genes %in% geneAnnot$name)){
        parallel::mclapply(mc.cores = nCores, genes, function(iGene){
            tmp2 <- geneAnnot[which(geneAnnot$name == iGene)]
            tmp3 <- c(GenomicAlignments::seqnames(tmp2), GenomicRanges::strand(tmp2)) %>%
                as.character() %>% c(iGene)
            exonBlock <- exonBlockGen(iGene, geneAnnot)
            tmpDT <- rep(tmp3, length(exonBlock)) %>%
                matrix(ncol = 3, byrow = T) %>% cbind(exonBlock) %>%
                cbind(seq(1, length(exonBlock)))
            tmpDT <- tmpDT[, c(1, 4, 2, 3, 5)] %>% data.table::data.table() %>%
                magrittr::set_colnames(c("chr", "gencoor", "strand",
                                         "gene", "txcoor"))
            tmpDT$chr <- as.factor(tmpDT$chr)
            tmpDT$gencoor <- as.integer(tmpDT$gencoor)
            tmpDT$strand <- as.factor(tmpDT$strand)
            tmpDT$gene <- as.factor(tmpDT$gene)
            tmpDT$txcoor <- as.integer(tmpDT$txcoor)
            if(!is.null(genome)){
                tmpDT <- tx_add_refSeqDT(tmpDT, genome, geneAnnot)
            }
            return(tmpDT)
        }) %>% magrittr::set_names(genes)
    }else{
        stop("Not all gene names are in geneAnnot")
    }
}

# Manipulating DTs #############################################################

# Removing UTR portions of DT, assuming first bases are
hlp_remove_UTR <- function(x, cut_5p = 0, cut_3p = 0){
    if(!all(diff(x$txcoor) == 1)){
        stop("Transcript coordinates 'txcoors' column is not continuous or has ",
             "gaps for gene: ", unique(x$gene))
    }
    tmp <- x[x$txcoor %in% (utils::tail(x$txcoor, -cut_5p) %>% utils::head(-cut_3p)),]
    tmp$txcoor <- tmp$txcoor - cut_5p
    return(tmp)
}

# Removing column if present
hlp_removeColumnIfPresent <- function(DT, colName){
    if(colName %in% colnames(DT)){
        DT[, colnames(DT) != colName, with = FALSE]
    }else{DT}
}

# Plotting helper funs #########################################################
txBrowser_colors <- list(
    "li_gray"   = "#d3d3da",
    "l_gray"   = "#b1b1be",
    "k_green" = "#00c201",
    "d_blue"    = "#0098fd",
    "h_yellow" = "#f3b018",
    "scarlet"  = "#fc2b06",
    "white" = "#ffffff",
    "r_black" = "#0b090b")
txBrowser2_colors <- list(
    "l_gray"   = "#b1b1be",
    "k_green" = "#00c201",
    "d_blue"    = "#0098fd",
    "h_yellow" = "#f3b018",
    "scarlet"  = "#fc2b06",
    "white" = "#ffffff",
    "r_black" = "#0b090b")

txBrowser_pal <- function(primary = "l_gray", other = "li_gray", direction = 1){
    # stopifnot(primary %in% names(txBrowser_colors))
    function(n) {
        if (n > 8) warning("txBrowser Color Palette only has 8 colors.")

        if (n == 2) {
            other <- if (!other %in% names(txBrowser_colors)) {
                other
            } else {
                txBrowser_colors[other]
            }
            color_list <- c(other, txBrowser_colors[primary])
        } else {
            color_list <- txBrowser_colors[1:n]
        }
        color_list <- unname(unlist(color_list))
        if (direction >= 0) color_list else rev(color_list)
    }
}
txBrowser_pal_2 <- function(direction = 1){
    function(n) {
        if (n > 7) warning("txBrowser Color Palette only has 8 colors.")
        color_list <- txBrowser2_colors[1:n]
        color_list <- unname(unlist(color_list))
        if (direction >= 0){color_list}else{rev(color_list)}
    }
}

# scale_fill with and without insert color
scale_fill_txBrowser <- function(primary = "l_gray", other = "li_gray", direction = 1) {
    ggplot2::discrete_scale("fill", "txBrowser",
                            txBrowser_pal(primary, other, direction))
}

scale_fill_txBrowser_2 <- function(direction = 1) {
    ggplot2::discrete_scale("fill", "txBrowser_2",
                            palette = txBrowser_pal_2(direction))
}

# Checking data objects and environment ########################################

# Check for data.table class, if DT is dataframe, convert to data.table
check_DT <- function(DT){
    if(!data.table::is.data.table(DT)){
        if(!is.data.frame(DT)){
            stop("DT must be either a data.table or a data.frame")
        }else{
            if(is.data.frame(DT)){
                DT <- data.table::data.table(DT)
            }
        }
    }
    return(DT)
}

# Check if OS is unix-like to allow multi-core operations
check_windows <- function(){Sys.info()[['sysname']] == "Windows"}

# Stop if nCores is greater than 1 and OS is windows
check_mc_windows <- function(nCores){
    if(nCores > 1 & check_windows()){
        stop("The multi-core capability of this function is ", "
             not available in Windows operating systems")
    }
}

# Check that DT object has $refSeq column
check_refSeq <- function(DT){
    if(!("refSeq" %in% names(DT))){
        stop("DT must contain the column 'refSeq' with the transcript ",
             "reference sequence")
    }
}

# Check that gene annotation and reads to be processed share chromosomes.
check_GA_reads_compatibility <- function(reads, geneAnnot){
    intChr <- intersect(as.character(unique(GenomeInfoDb::seqnames(geneAnnot))),
                        as.character(unique(GenomeInfoDb::seqnames(reads))))
    if(length(intChr) == 0){
        stop("reads and geneAnnot objects do not have chromosome names in common.")
    }
}

# Check that gene annotation and genome share chromosome names
check_GA_genome_chrCompat <- function(geneAnnot, genome){
    if(sum(unique(GenomeInfoDb::seqnames(geneAnnot)) %in% names(genome)) == 0){
        stop("No coinciding chromosome names between geneAnnot and genome")
    }
}

# Check integer argument
check_integer_arg <- function(x, argName){
    if(!is.numeric(x)){stop("Argument '", argName, "' must be an integer")}
    if((x - floor(x)) > 0){stop("Argument '", argName, "' must be an integer")}
}

# Check argument
check_integerGreaterThanZero_arg <- function(x, argName){
    check_integer_arg(x, argName)
    if(x <= 0){stop("Argument '", argName, "' must be greater than zero")}
}

# Check that resulting GRanges have the seq meta-column
check_GR_has_seq <- function(x, argName){
    if(class(x) == "CompressedGRangesList"){
        x <- unlist(x)
    }
    if(!("seq" %in% names(GenomicRanges::mcols(x)))){
        stop("'", argName, "' has no seq (sequences) meta-column.",
             "\nMake sure you set the argument 'withSeq' to TRUE in tx_reads()")
    }
}

check_BAM_has_seq <- function(x){
    pairedEnd <- switch(class(x), GAlignmentPairs = TRUE, GAlignments = FALSE)
    if(pairedEnd){
        if(!all(c("seq" %in% names(GenomicRanges::mcols(x@first)),
                  "seq" %in% names(GenomicRanges::mcols(x@last))))){
            stop("Input GAlignmentPairs object has no 'seq' metacolumn for either
                 first or last reads")
        }
    }else{
        if(!"seq" %in% names(GenomicRanges::mcols(x))){
            stop("Input GAlignments object has no 'seq' metacolumn for either ",
                 "first or last reads, required to work with nucleotide sequence. ",
                 "Load sequence information using the tx_load_bam() function by ",
                 "setting the 'loadSeq' argument to TRUE.")
        }
    }
}

#TODO: manage omitted ranges, print a warning if possible
# Check that end of range is >= to start -1, to make GenomicRanges from dataframes
check_DFforGRanges <- function(x){
    x[x$end >= (x$start -1),]
}

# Check that txDTs have same genes
check_sameGenesInDTL <- function(DTL){
    tmpL <- lapply(DTL, function(x) base::unique(x$gene))
    tmpA <- base::Reduce(x = tmpL, union)
    tmpB <- base::Reduce(x = tmpL, intersect)
    all(tmpA %in% tmpB)
}

# Annotate CDS start in tx_get_metageneAtCDS()
annot_CDSsta_DTL <- Vectorize(function(DT, CDS_start){
    DT$CDS_start[DT$gencoor == CDS_start[CDS_start$gene == DT$gene[1],]$gencoor] <- TRUE
    DT}, vectorize.args = "DT", SIMPLIFY = FALSE)

# Annotate CDS end. tx_get_metageneAtCDS()
annot_CDSend_DTL <- Vectorize(function(DT, CDS_end){
    DT$CDS_end[DT$gencoor == CDS_end[CDS_end$gene == DT$gene[1],]$gencoor] <- TRUE
    DT
    }, vectorize.args = "DT", SIMPLIFY = FALSE)


# Check that all genes in DT are in geneAnnot
check_GA_txDT_compat <- function(DT, geneAnnot){
    check_DThasCol(DT, "gene")
    if(!all(as.character(unique(DT$gene)) %in% geneAnnot$name)){
        stop("Not all genes in DT are contained in geneAnnot")
    }
}

# Check DT has a column with nam 'colName'
check_DThasCol <- function(DT, colName){
    DT <- check_DT(DT)
    if(!colName %in% colnames(DT)){stop("DT does not contain '", colName, "' column.")}
}

# "Hidden" functions (to decide if they'll be incorporated) ###################
#' Merge lists of data.tables
#'
#' @param DTL1 list List of data.table with gene names. As output of the
#' \code{\link{tx_coverageDT}}, \code{\link{tx_nucFreqDT}}, and
#' \code{\link{tx_covNucFreqDT}} functions.
#' @param DTL2 list List of data.table with gene names.
#' @param colsToAdd character. Numeric column(s) to be aggregated (added).
#' @param keepAll logical. Set to FALSE for just keeping data.tables which
#' name's are in both DTL.
#'
#' @return list
#'
#' @examples
hid_aggregate_DTlist <- function (DTL1, DTL2, colsToAdd, keepAll = TRUE){
    allNames <- union(names(DTL1), names(DTL2))
    namesInBoth <- intersect(names(DTL1), names(DTL2))
    namesOnly1 <- setdiff(names(DTL1), names(DTL2))
    namesOnly2 <- setdiff(names(DTL2), names(DTL1))
    tmpDTL <- lapply(namesInBoth, function(iGene){
        if(!identical(DTL1[[iGene]][, c("chr", "gencoor", "strand", "gene", "txcoor")],
                      DTL2[[iGene]][, c("chr", "gencoor", "strand", "gene", "txcoor")])){
            stop(paste("Coordinate data in", iGene, " data.tables are not compatible"))
        }else{
            tmp <- DTL1[[iGene]][, colsToAdd, with = FALSE] +
                DTL2[[iGene]][, colsToAdd, with = FALSE]
            cbind(DTL1[[iGene]][, c("chr", "gencoor", "strand", "gene", "txcoor"), with = FALSE], tmp)
        }
    }) %>% magrittr::set_names(namesInBoth)
    if(keepAll){
        c(tmpDTL, DTL1[namesOnly1], DTL2[namesOnly2])[allNames]
    }else{
        tmpDTL
    }
}

indexAlignmentsByGenomicRegion <- function(GAlignments, geneAnnot, overlapType = "within"){
    splitByChr <- split(geneAnnot, GenomicRanges::seqnames(geneAnnot))
    chrLen <- lapply(splitByChr, function(x) max(GenomicRanges::end(x))) %>% unlist()
    exonic <- unlist(exonGRanges(geneAnnot))
    allChr_GR <- rbind(data.frame(seqnames = names(chrLen),
                                  start = 1 ,
                                  end = chrLen,
                                  strand = "+"),
                       data.frame(seqnames = names(chrLen),
                                  start = 1 ,
                                  end = chrLen,
                                  strand = "-")) %>% plyranges::as_granges()
    notExon <- plyranges::setdiff_ranges_directed(allChr_GR, exonic)
    tmpGR <- GenomicRanges::findOverlaps(notExon, geneAnnot, type = overlapType)@from
    intronic <- notExon[tmpGR]
    intergen <- notExon[-tmpGR]
    if(class(GAlignments) == "GAlignmentPairs"){
        ovExonic <- union(GenomicRanges::findOverlaps(
            GAlignments@first, exonic, type = overlapType)@from,
            GenomicRanges::findOverlaps(
                GAlignments@last, GenomicRanges::invertStrand(exonic), type = overlapType)@from)
        ovIntron <- union(GenomicRanges::findOverlaps(
            GAlignments@first, intronic, type = overlapType)@from,
            GenomicRanges::findOverlaps(
                GAlignments@last, GenomicRanges::invertStrand(intronic), type = overlapType)@from)
        ovInterg <- union(GenomicRanges::findOverlaps(
            GAlignments@first, intergen, type = overlapType)@from,
            GenomicRanges::findOverlaps(
                GAlignments@last, GenomicRanges::invertStrand(intergen), type = overlapType)@from)
        ovRest <- setdiff(1:length(GAlignments), union(union(ovExonic, ovIntron), ovInterg))
    }else if(class(GAlignments) == "GAlignments"){
        ovExonic <- GenomicRanges::findOverlaps(GAlignments, exonic, type = overlapType)@from
        ovIntron <- GenomicRanges::findOverlaps(GAlignments, intronic, type = overlapType)@from
        ovInterg <- GenomicRanges::findOverlaps(GAlignments, intergen, type = overlapType)@from
        ovRest <- setdiff(1:length(GAlignments), union(union(ovExonic, ovIntron), ovInterg))
    }
    list(exonic = ovExonic,
         intronic = ovIntron,
         intergenic = ovInterg,
         rest = ovRest) %>% return()
}

# RowMeans by column groups
rowMeansColG <- function(DF, colGroups, na.rm = T){
    cG <- unique(colGroups)
    out <- sapply(cG, function(x){
        rowMeans(DF[,colGroups == x], na.rm = na.rm)
    }) %>% as.data.frame()
    colnames(out) <- cG
    rownames(out) <- rownames(DF)
    return(out)
}

#' Generate single-end FASTQ file
#'
#' Generate single-end FASTQ file form genome and gene annotation.
#' Distribution of reads is randomly selected following a negative binomial distribution with
#'
#' @param genome
#' @param geneAnnot
#' @param readLen
#' @param libSize
#' @param fileName
#' @param NB_r
#' @param NB_mu
#' @param nCores
#'
#' @return
#'
#' @examples
tx_generateSingleEndFASTQ <- function(genome, geneAnnot, readLen, libSize, fileName, NB_r = 5, NB_mu = 500, nCores){
    # filter transcripts by size
    txOME_seqs <- tx_get_transcriptSeqs(genome = genome, geneAnnot = geneAnnot, nCores = nCores)
    txOME_seqs <- txOME_seqs[BiocGenerics::width(txOME_seqs) >= readLen]
    # Random starts
    metaTX <- list(seqs = txOME_seqs, width = BiocGenerics::width(txOME_seqs))
    metaTX$TPM <- stats::rnbinom(n = length(txOME_seqs), size = NB_r, mu = NB_mu)
    metaTX$TPM <- round(metaTX$TPM * (libSize / sum(metaTX$TPM)))
    metaTX$range1 <- metaTX$width - readLen + 1
    metaTX$starts_1 <- lapply(seq_along(metaTX$TPM), function(i){
        sample(seq(1, metaTX$range1[i]), size = metaTX$TPM[i], replace = TRUE)
    })
    #Extract sequences
    metaTX$reads1 <- lapply(seq_along(metaTX$seqs), function(i){
        stringr::str_sub(metaTX$seqs[i], metaTX$starts_1[[i]],
                         metaTX$starts_1[[i]] + readLen - 1) %>%
            Biostrings::DNAStringSet()
    }) %>% do.call(what = "c")
    names(metaTX$reads1) <- paste0("R1_", 1:length(metaTX$reads1))
    qualStr <- paste(rep("H", readLen), collapse = "")
    # Writing FASTQ
    Biostrings::writeXStringSet(x = metaTX$reads1,
                                quali = Biostrings::BStringSet(rep(qualStr, each = length(metaTX$reads1))),
                                filepath = fileName, format = "fastq", compress = TRUE)
}


## usethis namespace: start
#' @importFrom lifecycle deprecated
## usethis namespace: end

