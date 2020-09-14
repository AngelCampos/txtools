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
    listLen <- lapply(tmpA, length) %>% unlist
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
    split_i <- split_i[sapply(split_i, length) %>%
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
    # Filtering reads extrictly inside exons
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
        tReads <- data.frame(start = match(GenomicRanges::start(iReads_r1[pass]), iExon),
                             end      = match(GenomicRanges::end(iReads_r2[pass]), iExon),
                             strand   = GenomicRanges::strand(iReads_r1[pass]),
                             seqnames = iGene) %>% plyranges::as_granges()
    }else if(iStrand == "-"){
        tReads <- data.frame(start = match(GenomicRanges::end(iReads_r1[pass]), iExon),
                             end      = match(GenomicRanges::start(iReads_r2[pass]), iExon),
                             strand   = GenomicRanges::strand(iReads_r1[pass]),
                             seqnames = iGene) %>% plyranges::as_granges()
    }
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
                                                        S4Vectors::mcols(tReads)$cigar_r1) %>%
            stringr::str_remove_all(pattern = "\\.")
        tReads$seq2 <- GenomicAlignments::sequenceLayer(S4Vectors::mcols(tReads)$seq_r2,
                                                        S4Vectors::mcols(tReads)$cigar_r2) %>%
            stringr::str_remove_all(pattern = "\\.")
    }else if(iStrand == "-"){
        tReads$seq1 <- GenomicAlignments::sequenceLayer(S4Vectors::mcols(tReads)$seq_r1,
                                                        S4Vectors::mcols(tReads)$cigar_r1) %>%
            Biostrings::reverseComplement() %>% stringr::str_remove_all(pattern = "\\.")
        tReads$seq2 <- GenomicAlignments::sequenceLayer(S4Vectors::mcols(tReads)$seq_r2,
                                                        S4Vectors::mcols(tReads)$cigar_r2) %>%
            Biostrings::reverseComplement() %>% stringr::str_remove_all(pattern = "\\.")
    }
    # Trim overflowing reads
    i <- which(nchar(tReads$seq1) > GenomicAlignments::width(tReads))
    tReads[i]$seq1 <- stringr::str_sub(tReads[i]$seq1, start = 1, GenomicAlignments::width(tReads[i]))
    i <- which(nchar(tReads$seq2) > GenomicAlignments::width(tReads))
    tReads[i]$seq2 <- stringr::str_sub(tReads[i]$seq2, start = - GenomicAlignments::width(tReads[i]))
    # Calculate overlap
    tReads$diffSeq <- nchar(tReads$seq1) + nchar(tReads$seq2) - GenomicAlignments::width(tReads)
    tReads$oL <- tReads$diffSeq %>% magrittr::is_greater_than(0)
    tReads$mergedSeq <- "." # Place holder
    # No overlap reads
    i <- which(!(tReads$oL))
    if(length(i) > 0){
        gapLen <- GenomicAlignments::width(tReads[i]) - nchar(tReads[i]$seq1) - nchar(tReads[i]$seq2)
        insSeq <- stringr::str_dup(string = ".", gapLen)
        tReads[i]$mergedSeq <- paste(tReads[i]$seq1, insSeq,
                                     tReads[i]$seq2, sep = "") %>% Biostrings::DNAStringSet()
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
                    Biostrings::consensusString(ambiguityMap = txtools::IUPAC_CODE_MAP_extended, threshold = 0.2)
            }
        }
        tReads[i]$mergedSeq <- paste0(stringr::str_sub(tReads[i]$seq1,
                                                       start = 1,
                                                       end = nchar(tReads[i]$seq1) - tReads[i]$diffSeq),
                                      ovSeq, stringr::str_sub(tReads[i]$seq2,
                                                              start = tReads[i]$diffSeq + 1,
                                                              end = nchar(tReads[i]$seq2)))
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
    # Both ends of both reads fall into the gene model
    stEndTable <- as.matrix(data.frame(r1_S = GenomicRanges::start(iReads_r1),
                                       r1_E = GenomicRanges::end(iReads_r1)))
    pass <- (stEndTable %in% iExon) %>% matrix(ncol = 2, byrow = F) %>%
        rowSums %>% magrittr::equals(2) %>% which()
    if(length(pass) < minReads){return(GenomicRanges::GRanges())} # No reads Return empty GA
    # Reads cover consecutive exons
    tmp <- GenomicRanges::findOverlaps(iReads_r1[pass], allExons[[iGene]])
    passPos <- split(tmp@to, tmp@from) %>% lapply(diff) %>%
        lapply(function(x) all(x == 1)) %>% unlist %>% which
    pass <- pass[passPos]
    if(length(pass) < minReads){return(GenomicRanges::GRanges())} # No reads Return empty GA
    # Boundaries of merged reads
    if(iStrand == "+"){
        tReads <- data.frame(start = match(GenomicRanges::start(iReads_r1[pass]), iExon),
                             end      = match(GenomicRanges::end(iReads_r1[pass]), iExon),
                             strand   = GenomicRanges::strand(iReads_r1[pass]),
                             seqnames = iGene) %>% plyranges::as_granges()
    }else if(iStrand == "-"){
        tReads <- data.frame(start = match(GenomicRanges::end(iReads_r1[pass]), iExon),
                             end      = match(GenomicRanges::start(iReads_r1[pass]), iExon),
                             strand   = GenomicRanges::strand(iReads_r1[pass]),
                             seqnames = iGene) %>% plyranges::as_granges()
    }
    names(tReads) <- names(reads)[selReadsbyPair][pass]
    GenomeInfoDb::seqlengths(tReads) <- length(iExon)
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
                                                        S4Vectors::mcols(tReads)$cigar_r1) %>%
            stringr::str_remove_all(pattern = "\\.")
    }else if(iStrand == "-"){
        tReads$seq1 <- GenomicAlignments::sequenceLayer(S4Vectors::mcols(tReads)$seq_r1,
                                                        S4Vectors::mcols(tReads)$cigar_r1) %>%
            Biostrings::reverseComplement() %>%
            stringr::str_remove_all(pattern = "\\.")
    }
    # Trim overflowing reads
    i <- which(nchar(tReads$seq1) > GenomicAlignments::width(tReads))
    tReads[i]$seq1 <- stringr::str_sub(tReads[i]$seq1, start = 1,
                                       GenomicAlignments::width(tReads[i]))
    # Calculate overlap
    tReads$diffSeq <- nchar(tReads$seq1) - GenomicAlignments::width(tReads)
    tReads$oL <- tReads$diffSeq %>% magrittr::is_greater_than(0)
    tReads$mergedSeq <- "." # Place holder
    # No overlap reads
    i <- which(!(tReads$oL))
    if(length(i) > 0){
        gapLen <- GenomicAlignments::width(tReads[i]) - nchar(tReads[i]$seq1)
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

# Calculate coverage for all gene models
tx_coverage <- function(x){
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

# Vectorized intersect
vIntersect <- Vectorize(intersect, c("x", "y"), SIMPLIFY = F)

# Pipe
`%>%` <- magrittr::`%>%`

# Object size
oSize <- function(x){
    print(utils::object.size(x), units = "auto")
}

# Generating DTs ###############################################################
# Generates trasncriptomic coordinates table from a list of genes
hlpr_genCoorTabGenes <- function(genes, geneAnnot, fastaGenome = NULL, nCores = 1){
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
            if(!is.null(fastaGenome)){
                tmpDT <- tx_add_refSeqDT(tmpDT, fastaGenome, geneAnnot)
            }
            return(tmpDT)
        }) %>% magrittr::set_names(genes)
    }else{
        stop("Not all gene names are in geneAnnot")
    }
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

# Other ####

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
stop_mc_windows <- function(nCores){
    if(nCores > 1 & check_windows()){
        stop("The multi-core capability of this function is ", "
             not available in Windows operating systems.\n")
    }
}

# Check that DT object has $refSeq column
check_refSeq <- function(DT){
    if(!("refSeq" %in% names(DT))){
        stop("DT must contain the column 'refSeq' with the transcript ",
             "reference sequence")
    }
}
