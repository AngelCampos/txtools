#' Read paired end bam file by yield size
#'
#' @param file character. The path to the file to read
#' @param yieldSize numeric. Size of chunck of reads processed at a time
#' @param scanFlag integer. Flag used to filter reads check: ?scanBamFlag()
#' @param loadSeq logical. Set to TRUE for loading the sequences contained
#' in the BAM file
#' @param recoverDumpedAligns logical. If set to TRUE ambiguous alignments
#' will be retrieved and output along correctly paired reads. More info:
#' ?GenomicAlignments::readGAlignments
#' @param verbose logical. Set to FALSE to show less information
#'
#' @return
#' @export
#'
#' @examples
load_pairedEnd_bam <- function(file,
                               yieldSize = "auto",
                               scanFlag = "default",
                               loadSeq = FALSE,
                               recoverDumpedAligns = F,
                               verbose = TRUE){
    if(verbose){cat("Reading number of records in file \n")}
    #Remove optical duplicates, low quality reads, low quality reads, non-paired
    if(scanFlag == "default"){
        scanFlag <- Rsamtools::scanBamFlag(isDuplicate = FALSE,
                                           isNotPassingQualityControls = FALSE,
                                           isPaired = TRUE)
    }
    # CountBam records (each paired end alignments equals two records)
    bC <- Rsamtools::countBam(file, )
    if(bC$records == 0){stop("BAM file is empty \n")}
    if(verbose){cat(bC$records, "number of BAM records \n")}
    if(yieldSize == "auto"){
        yieldSize <- ceiling(bC$records/25)
    }
    if(verbose){cat("Loading BAM file \n")}
    pb <- txtProgressBar(style = 3)
    BAMFILE <- Rsamtools::BamFile(file, yieldSize = yieldSize)
    readCycles <- 1:ceiling(bC$records/(yieldSize*2))
    pbJumps <- seq(0,1, by = min(1/(length(readCycles)-1), 1))
    if(length(readCycles) > 200){stop("Too low yieldSize")}
    #  Opens and reads BAM files in chunks of length = yieldSize
    open(BAMFILE)
    reads <- lapply(pbJumps, function(i){
        if(verbose){setTxtProgressBar(pb, i)}
        if(loadSeq){
            GenomicAlignments::readGAlignmentPairs(BAMFILE, use.names = TRUE,
                                                   param = Rsamtools::ScanBamParam(flag = scanFlag, what = "seq"))
        }else if(loadSeq == F){
            GenomicAlignments::readGAlignmentPairs(BAMFILE, use.names = TRUE,
                                                   param = Rsamtools::ScanBamParam(flag = scanFlag))
        }else(stop("loadSeq parameter must be logical either FALSE or TRUE"))
    }) %>% do.call(what = c)
    close(BAMFILE)
    close(pb)
    bamData <- list(GAligns = reads, dumpedAmbigPairs = GenomicAlignments::getDumpedAlignments())
    if(verbose){cat(" \n")}
    if(verbose){cat(length(bamData$GAligns), "reads succesfully loaded \n")}
    cat("Dumped reads due to ambiguous pairs:", length(bamData$dumpedAmbigPairs), "\n")
    if(recoverDumpedAligns == F){
        return(reads)
    }else{
        return(bamData)
    }
}

#
#' Representing gene models as GRanges from imported BEDFILE
#'
#' @param geneAnnot_GR
#'
#' @return
#' @export
#'
#' @examples
exonGRanges <- function(geneAnnot_GR){
    iChr <- seqnames(geneAnnot_GR) %>% as.character()
    iStart <- BiocGenerics::start(geneAnnot_GR)
    iEnd <- BiocGenerics::end(geneAnnot_GR)
    iStrand <- strand(geneAnnot_GR)
    tmpA <- BiocGenerics::start(geneAnnot_GR$blocks) %>% subtract(1) %>% add(iStart)
    tmpB <- width(geneAnnot_GR$blocks) %>% subtract(1) %>% add(tmpA)
    listLen <- lapply(tmpA, length) %>% unlist
    group <- rep(1:length(geneAnnot_GR), times = listLen)
    data.frame(rep(iChr, times = listLen),
               unlist(tmpA),
               unlist(tmpB),
               rep(iStrand, times = listLen)) %>%
        set_colnames(c("seqnames", "start", "end", "strand")) %>%
        as_granges() %>%
        split(group) %>%
        GRangesList() %>% set_names(geneAnnot_GR$name)
}

# Generate exon coordinates block
exonBlockGen <- function(iGene, geneAnnot_GR){
    iStart <- BiocGenerics::start(geneAnnot_GR[geneAnnot_GR$name == iGene])
    iEnd <- BiocGenerics::end(geneAnnot_GR[geneAnnot_GR$name == iGene])
    iStrand <- strand(geneAnnot_GR[geneAnnot_GR$name == iGene]) %>% as.character()
    tmpA <- unlist(BiocGenerics::start(geneAnnot_GR[geneAnnot_GR$name == iGene]$blocks)) -1
    tmpB <- unlist(width(geneAnnot_GR[geneAnnot_GR$name == iGene]$blocks))
    iBlocks <- sapply(1:length(tmpA), function(i){
        c(tmpA[i], (tmpA[i] + tmpB[i] -1))
    }) %>% t %>% add(iStart)
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

# Merge paired-end reads and assign them to gene models (Multicore)
tx_PEreads_mc <- function(reads, bedR, nCores, overlapType = "within",
                          minReads = 50, withSeq = F, verbose = T){
    if(class(reads) != "GAlignmentPairs"){
        stop("reads argument should be of class GAlignmentPairs \n")
    }
    if(class(bedR) != "GRanges"){
        stop("bedR argument should be of class GRanges \n")
    }
    if(.Platform$OS.type != "unix"){
        stop("This functions is only available for UNIX operative systems \n")
    }
    if(verbose){
        cat("Processing", length(reads), "paired-end reads, using", length(bedR), "gene models \n")
    }
    split_i <- hlpr_splitReadsByGenes(reads, bedR, overlapType,
                                      minReads)
    bedR <- bedR[which(bedR$name %in% names(split_i))]
    if(length(bedR) > 0){
        if(verbose){
            cat(length(unique(unlist(split_i))), "paired-end reads overlap",
                length(bedR), "gene models \n")
            cat("Filtering reads by gene model... \n")
            if(withSeq){
                cat("Processing sequences. This may take several minutes... \n")
            }
        }
    }else{stop("No genes with overlapped paired-end reads \n")}
    allExons <- exonGRanges(bedR) # All exons in gene models
    OUT <- mclapply(mc.cores = nCores, bedR$name, function(iGene){
        hlpr_ReadsInGene(reads = reads,
                         iGene = iGene,
                         bedR = bedR,
                         split_i = split_i,
                         allExons = allExons,
                         withSeq = withSeq,
                         minReads = minReads)
    })
    names(OUT) <- bedR$name
    OUT <- OUT[lapply(OUT, length) %>% unlist %>% is_greater_than(minReads)] %>%
        GenomicRangesList()
    if(verbose){
        cat("Output contains:", lapply(OUT, names) %>% unlist %>% unique %>% length,
            "unique reads in", length(OUT), "gene models \n")
    }
    return(OUT)
}

#' Merge paired-end reads and assign them to gene models
#'
#' @param reads
#' @param bedR
#' @param overlapType
#' @param minReads
#' @param withSeq
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
tx_PEreads <- function(reads, bedR, overlapType = "within", minReads = 50, withSeq = F, verbose = T){
    if(class(reads) != "GAlignmentPairs"){
        stop("reads argument should be of class GAlignmentPairs \n")
    }
    if(class(bedR) != "GRanges"){
        stop("bedR argument should be of class GRanges \n")
    }
    if(.Platform$OS.type != "unix"){
        stop("This functions is only available for UNIX operative systems \n")
    }
    if(verbose){
        cat("Processing", length(reads), "reads, using", length(bedR), "gene models \n")
    }
    split_i <- hlpr_splitReadsByGenes(reads, bedR, overlapType, minReads)
    bedR <- bedR[which(bedR$name %in% names(split_i))]
    if(length(bedR) > 0){
        if(verbose){
            cat(length(unique(unlist(split_i))), "paired-end reads overlap",
                length(bedR), "gene models \n")
            cat("Filtering reads by gene model... \n")
            if(withSeq){
                cat("Processing sequences. This may take several minutes... \n")
            }
        }
    }else{stop("No genes with overlapped paired-end reads \n")}
    allExons <- exonGRanges(bedR) # All exons in gene models
    OUT <- v_hlpr_ReadsInGene(reads = reads,
                              iGene = bedR$name,
                              bedR = bedR,
                              split_i = split_i,
                              allExons = allExons,
                              withSeq = withSeq,
                              minReads = minReads)
    names(OUT) <- bedR$name
    OUT <- OUT[lapply(OUT, length) %>% unlist %>% is_greater_than(minReads)] %>%
        GenomicRangesList()
    if(verbose){
        cat("Output contains:", lapply(OUT, names) %>% unlist %>% unique %>% length,
            "unique reads in", length(OUT), "gene models \n")
    }
    return(OUT)
}

# Reads overlapping gene models
hlpr_splitReadsByGenes <- function(reads, bedR, overlapType, minReads){
    #IgnoreStrand feature missing
    allOver_1 <- GenomicRanges::findOverlaps(reads@first, bedR, type = overlapType)
    allOver_2 <- GenomicRanges::findOverlaps(GenomicAlignments::invertStrand(reads@last), bedR,type = overlapType)
    split_1 <- split(allOver_1@from, allOver_1@to)
    split_2 <- split(allOver_2@from, allOver_2@to)
    names(split_1) <- bedR[names(split_1) %>% as.numeric()]$name
    names(split_2) <- bedR[names(split_2) %>% as.numeric()]$name
    inBoth <- intersect(names(split_1), names(split_2))
    split_1 <- split_1[inBoth]
    split_2 <- split_2[inBoth]
    split_i <- vIntersect(split_1, split_2)
    names(split_i) <- names(split_1)
    split_i <- split_i[sapply(split_i, length) %>%
                           magrittr::is_weakly_greater_than(minReads) %>% which]
    return(split_i)
}

hlpr_ReadsInGene <- function(reads, iGene, bedR, split_i, allExons, minReads, withSeq){
    iStrand <- bedR[which(bedR$name == iGene)] %>% strand %>% as.character()
    iExon <- exonBlockGen(iGene = iGene, geneAnnot_GR = bedR)
    selReadsbyPair <- split_i[[iGene]]
    # Selecting paired reads to merge
    iReads_r1 <- reads@first[selReadsbyPair]
    iReads_r2 <- reads@last[selReadsbyPair]
    # Filtering reads extrictly inside exons
    # Both ends of both reads fall into the gene model
    stEndTable <- as.matrix(data.frame(r1_S = BiocGenerics::start(iReads_r1),
                                       r1_E = BiocGenerics::end(iReads_r1),
                                       r2_S = BiocGenerics::start(iReads_r2),
                                       r2_E = BiocGenerics::end(iReads_r2)))
    pass <- (stEndTable %in% iExon) %>% matrix(ncol = 4, byrow = F) %>%
        rowSums %>% equals(4) %>% which
    if(length(pass) < minReads){return(GRanges())} # No reads Return empty GA
    # Reads cover consecutive exons
    tmp <- GenomicRanges::findOverlaps(iReads_r1[pass], allExons[[iGene]])
    passPos <- split(tmp@to, tmp@from) %>% lapply(diff) %>%
        lapply(function(x) all(x == 1)) %>% unlist %>% which
    tmp <- GenomicRanges::findOverlaps(GenomicAlignments::invertStrand(iReads_r2[pass]),
                                       allExons[[iGene]])
    passNeg <- split(tmp@to, tmp@from) %>% lapply(diff) %>%
        lapply(function(x) all(x == 1)) %>% unlist %>% which
    pass <- pass[as.numeric(intersect(passNeg, passPos))]
    if(length(pass) < minReads){return(GRanges())} # No reads Return empty GA
    # Boundaries of merged reads
    if(iStrand == "+"){
        tReads <- data.frame(start = match(BiocGenerics::start(iReads_r1[pass]), iExon),
                             end      = match(BiocGenerics::end(iReads_r2[pass]), iExon),
                             strand   = strand(iReads_r1[pass]),
                             seqnames = iGene) %>% as_granges
    }else if(iStrand == "-"){
        tReads <- data.frame(start = match(BiocGenerics::end(iReads_r1[pass]), iExon),
                             end      = match(BiocGenerics::start(iReads_r2[pass]), iExon),
                             strand   = strand(iReads_r1[pass]),
                             seqnames = iGene) %>% as_granges
    }
    names(tReads) <- names(reads)[selReadsbyPair][pass]
    seqlengths(tReads) <- length(iExon)
    # Result if no sequence is input or required
    if(withSeq == F){
        return(tReads)
    }
    # Merging sequences
    tReads$start_r1 <- BiocGenerics::start(iReads_r1[pass])
    tReads$end_r1 <- BiocGenerics::end(iReads_r1[pass])
    tReads$start_r2 <- BiocGenerics::start(iReads_r2[pass])
    tReads$end_r2 <- BiocGenerics::end(iReads_r2[pass])
    tReads$strand_r1 <- strand(iReads_r1[pass])
    tReads$cigar_r1 <- cigar(iReads_r1[pass])
    tReads$cigar_r2 <- cigar(iReads_r2[pass])
    tReads$seq_r1 <- mcols(iReads_r1[pass])$seq
    tReads$seq_r2 <- mcols(iReads_r2[pass])$seq
    # Constructing the merged read sequence
    if(iStrand == "+"){
        tReads$seq1 <- sequenceLayer(mcols(tReads)$seq_r1,
                                     mcols(tReads)$cigar_r1) %>%
            str_remove_all(pattern = "\\.")
        tReads$seq2 <- sequenceLayer(mcols(tReads)$seq_r2,
                                     mcols(tReads)$cigar_r2) %>%
            str_remove_all(pattern = "\\.")
    }else if(iStrand == "-"){
        tReads$seq1 <- sequenceLayer(mcols(tReads)$seq_r1,
                                     mcols(tReads)$cigar_r1) %>%
            reverseComplement() %>% str_remove_all(pattern = "\\.")
        tReads$seq2 <- sequenceLayer(mcols(tReads)$seq_r2,
                                     mcols(tReads)$cigar_r2) %>%
            reverseComplement() %>% str_remove_all(pattern = "\\.")
    }
    # Trim overflowing reads
    i <- which(nchar(tReads$seq1) > width(tReads))
    tReads[i]$seq1 <- str_sub(tReads[i]$seq1, start = 1, width(tReads[i]))
    i <- which(nchar(tReads$seq2) > width(tReads))
    tReads[i]$seq2 <- str_sub(tReads[i]$seq2, start = -width(tReads[i]))
    # Calculate overlap
    tReads$diffSeq <- nchar(tReads$seq1) + nchar(tReads$seq2) - width(tReads)
    tReads$oL <- tReads$diffSeq %>% is_greater_than(0)
    tReads$mergedSeq <- "." # Place holder
    # No overlap reads
    i <- which(!(tReads$oL))
    if(length(i) > 0){
        gapLen <- width(tReads[i]) - nchar(tReads[i]$seq1) - nchar(tReads[i]$seq2)
        insSeq <- str_dup(string = ".", gapLen)
        tReads[i]$mergedSeq <- paste(tReads[i]$seq1, insSeq,
                                     tReads[i]$seq2, sep = "") %>%
            DNAStringSet()
    }
    # Overlapped reads
    i <- which(tReads$oL)
    if(length(i) > 0){
        tmp1 <- str_sub(tReads[i]$seq1, start = -tReads[i]$diffSeq)
        tmp2 <- str_sub(tReads[i]$seq2, start = 1, end = tReads[i]$diffSeq)
        ovSeq <- rep(".", length(i))
        for(j in 1:length(i)){
            if(tmp1[j] == tmp2[j]){
                ovSeq[j] <- tmp1[j]
            }else{
                ovSeq[j] <- DNAStringSet(c(tmp1[j], tmp2[j])) %>%
                    consensusMatrix %>%
                    consensusString(ambiguityMap = IUPAC_CODE_MAP_extended, threshold = 0.2)
            }
        }
        tReads[i]$mergedSeq <- paste0(str_sub(tReads[i]$seq1, start = 1, end = nchar(tReads[i]$seq1) - tReads[i]$diffSeq),
                                      ovSeq,
                                      str_sub(tReads[i]$seq2, start = tReads[i]$diffSeq + 1, end = nchar(tReads[i]$seq2)))
    }
    # Final reads
    fReads <- tReads
    mcols(fReads) <- NULL
    fReads$seq <- tReads$mergedSeq
    return(fReads)
}

# Vectorized version of helper
v_hlpr_ReadsInGene <- Vectorize(hlpr_ReadsInGene, "iGene")

# Read counts by gene model
tx_counts <- function(x){
    suppressWarnings(unlist(x)) %>% seqnames %>% table
}

# Calculate coverage for all gene models
tx_coverage <- function(x){
    suppressWarnings(unlist(x)) %>% GenomicRanges::coverage()
}

# Calculate coverage table (cov, ends, starts) for all gene models
tx_covTab <- function(x){
    if(class(x) != "SimpleGenomicRangesList"){
        stop("x must be of class SimpleGenomicRangesList")
    }
    if(names(x) %>% duplicated() %>% sum() %>% magrittr::is_greater_than(0)){
        stop("List contains duplicated gene models")
    }
    cov <- tx_coverage(x) %>% lapply(as.vector)
    sta <- BiocGenerics::start(x) %>% lapply(table) %>% set_names(names(x))
    end <- BiocGenerics::end(x) %>% lapply(table) %>% set_names(names(x))
    len <- lapply(x, seqlengths) %>% unlist
    OUT <- lapply(1:length(x), function(i){
        tmpMat <- matrix(0, nrow = len[i], ncol = 3) %>% data.frame() %>% set_rownames(1:len[i]) %>%
            set_colnames(c("cov", "start_5p", "end_3p"))
        tmpMat[names(sta[[i]]), "start_5p"] <- sta[[i]]
        tmpMat[names(end[[i]]), "end_3p"] <- end[[i]]
        tmpMat[, "cov"] <- cov[[i]] %>% as.numeric()
        # tmpMat[,"start_5p"] <- Rle(tmpMat[,"start_5p"])
        # tmpMat[,"end_3p"] <- Rle(tmpMat[,"end_3p"])
        tmpMat %>% data.table()
    }) %>% set_names(names(x))
    return(OUT)
}

# Calculate coverage table (cov, ends, starts) for all gene models (Multicore)
tx_covTab_mc <- function(x, nCores){
    if(class(x) != "SimpleGenomicRangesList"){
        stop("x must be of class SimpleGenomicRangesList")
    }
    if(names(x) %>% duplicated() %>% sum() %>% magrittr::is_greater_than(0)){
        stop("List contains duplicated gene models")
    }
    cov <- tx_coverage(x) %>% mclapply(as.vector, mc.cores = nCores)
    sta <- BiocGenerics::start(x) %>% lapply(table) %>% set_names(names(x))
    end <- BiocGenerics::end(x) %>% lapply(table) %>% set_names(names(x))
    len <- lapply(cov, length) %>% unlist
    OUT <- mclapply(mc.cores = nCores, 1:length(x), function(i){
        tmpMat <- matrix(0, nrow = len[i], ncol = 3) %>% data.frame() %>% set_rownames(1:len[i]) %>%
            set_colnames(c("cov", "start_5p", "end_3p"))
        tmpMat[names(sta[[i]]), "start_5p"] <- sta[[i]]
        tmpMat[names(end[[i]]), "end_3p"] <- end[[i]]
        tmpMat[, "cov"] <- cov[[i]] %>% as.numeric()
        # tmpMat[,"start_5p"] <- Rle(tmpMat[,"start_5p"])
        # tmpMat[,"end_3p"] <- Rle(tmpMat[,"end_3p"])
        tmpMat %>% data.table()
    }) %>% set_names(names(x))
    return(OUT)
}

# Filter by width
tx_filter_max_width <- function(x, thr){
    tmp <- GenomicAlignments::width(x) %>% magrittr::is_weakly_less_than(thr)
    lapply(seq(1, length(x)), function(i){
        if(all(tmp[[i]])){
            x[[i]]
        }else{
            x[[i]][tmp[[i]]]
        }
    }) %>% GenomicRangesList() %>% set_names(names(x))
}

# Calculate nucleotide frequency pileup for all gene models

tx_nucFreqTab <- function(x, simplify_IUPAC = "not"){
    OUT <- lapply(seq(1, length(x)), function(i){
        y <- Biostrings::consensusMatrix(x[[i]]$seq, shift = BiocGenerics::start(x[[i]]) -1,
                                         width = seqlengths(x[[i]]))
        if(simplify_IUPAC == "not"){
            hlp_addMissingNucs(y) %>% t %>% data.table()
        }else if(simplify_IUPAC == "splitHalf"){
            hlp_splitNucsHalf(y) %>% t %>% data.table()
        }else if(simplify_IUPAC == "splitForceInt"){
            hlp_splitNucsForceInt(y) %>% t %>% data.table()
        }
    }) %>% set_names(names(x))
}

hlp_addMissingNucs <- function(x){
    misNucs <- IUPAC_code_2nucs[which(!(IUPAC_code_2nucs %in% rownames(x)))]
    if(length(misNucs) > 0){
        matrix(0, nrow = length(misNucs), ncol = ncol(x)) %>%
            set_rownames(misNucs) %>% rbind(x) %>% .[IUPAC_code_2nucs,]
    }else{
        x
    }
}

hlp_splitNucsHalf <- function(x){
    misNucs <- IUPAC_code_simpl[which(!(IUPAC_code_simpl %in% rownames(x)))]
    altNucs <- intersect(IUPAC_code_2nucs[5:10], rownames(x))
    x <- hlp_addMissingNucs(x)
    for(i in altNucs){
        resNuc <- IUPAC_CODE_MAP[i] %>% str_split("") %>% unlist
        x[resNuc[1],] <- x[resNuc[1],] + (x[i, ] / 2)
        x[resNuc[2],] <- x[resNuc[2],] + (x[i, ] / 2)
    }
    x[IUPAC_code_simpl,]
}

hlp_splitNucsForceInt <- function(x){
    misNucs <- IUPAC_code_simpl[which(!(IUPAC_code_simpl %in% rownames(x)))]
    altNucs <- intersect(IUPAC_code_2nucs[5:10], rownames(x))
    x <- hlp_addMissingNucs(x)
    for(i in altNucs){
        if(all(x[i,] %% 2 == 0)){
            resNuc <- IUPAC_CODE_MAP[i] %>% str_split("") %>% unlist
            x[resNuc[1],] <- x[resNuc[1],] + (x[i,] / 2)
            x[resNuc[2],] <- x[resNuc[2],] + (x[i,] / 2)
        }else{
            resNuc <- IUPAC_CODE_MAP[i] %>% str_split("") %>% unlist
            x[resNuc[1],] <- x[resNuc[1],] + (x[i,] %>% divide_by(2) %>% floor)
            x[resNuc[2],] <- x[resNuc[2],] + (x[i,] %>% divide_by(2) %>% floor)
            x["N", ] <- x["N",] + (x[i,] - (x[i,] %>% divide_by(2) %>% floor %>% multiply_by(2)))
        }
    }
    x[IUPAC_code_simpl,]
}

# Table with genomic and transcriptomic coordinates
tx_genCoorTab <- function(x, geneAnnot_GR){
    if(all(names(x) %in% geneAnnot_GR$name)){
        lapply(names(x), function(iGene){
            tmp2 <- geneAnnot_GR[which(geneAnnot_GR$name == iGene)]
            tmp3 <- c(seqnames(tmp2), strand(tmp2)) %>% as.character() %>% c(iGene)
            rep(tmp3, seqlengths(x[[iGene]])) %>% matrix(ncol = 3, byrow = T) %>%
                cbind(exonBlockGen(iGene, geneAnnot_GR)) %>%
                cbind(seq(1, seqlengths(x[[iGene]]))) %>% .[,c(1,4,2,3,5)] %>%
                data.table() %>% set_colnames(c("chr", "gencoor", "strand", "gene", "txcoor"))
        }) %>% set_names(names(x))
    }else{
        stop("Names of x are not contained in geneAnnot_GR$name")
    }
}

# Helper bind three results tables
hlp_cbind3Tabs <- function(gencoorT, tab1, tab2){
    if(all(names(gencoorT) == names(tab1) &
           names(tab2) == names(gencoorT))){
        lapply(seq(1, length(gencoorT)), function(i){
            cbind(gencoorT[[i]], tab1[[i]], tab2[[i]])
        }) %>% set_names(names(gencoorT))
    }
}
# Helper bind two results tables
hlp_cbind2Tabs <- function(gencoorT, tab1){
    if(all(names(gencoorT) == names(tab1))){
        lapply(seq(1, length(gencoorT)), function(i){
            cbind(gencoorT[[i]], tab1[[i]])
        }) %>% set_names(names(gencoorT))
    }
}

# Vectorized intersect
vIntersect <- Vectorize(intersect, c("x", "y"), SIMPLIFY = F)

