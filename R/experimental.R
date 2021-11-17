tx_plot_metaGeneByBins <- function(DT, colName, nBins = 100, FUN = "mean", minTxLength = 300,
                                   minReadsPerGene = 100, smooth = TRUE, na.rm = TRUE){
    if(nBins >= minTxLength){stop("Number of bins most be smaller than minimum",
                                  "transcript length.")}
    if(!(colName %in% colnames(DT))){stop("colName is not a column in DT.")}
    tapply(txDTn[["start_5p"]], txDT$gene, sum)
    nStarts <- DT[, sum(start_5p), by = gene][order(V1, decreasing = T),]
    DT <- DT[gene %in% nStarts$gene[nStarts$V1 >= minReadsPerGene]]
    geneLens <- tx_get_geneLengths(DT)
    DT <- DT[gene %in% names(geneLens)[geneLens > minTxLength]]
    GENES <- as.character(unique(DT$gene))
    if(length(GENES) < 1){stop("No genes after filtering parameters")}
    meanCovBinned <- mclapply(seq(GENES), function(i){
        iGene <- GENES[i]
        tmpDT <- DT[gene == iGene,]
        tmpDT$group <- tmpDT$txcoor %>% cut_number(n = nBins)
        out <- tapply(tmpDT[[colName]], tmpDT$group, FUN, na.rm = na.rm)
        out[is.nan(out)] <- NA
        out
    }) %>% do.call(what = rbind) %>% apply(MARGIN = 2, FUN = FUN, na.rm = na.rm) %>% set_names(NULL)
    DF <- data.frame(bins = seq(meanCovBinned) %>% as.numeric,
                     score = meanCovBinned)
    plotTitle <- paste("METAGENE BY BINS -", FUN, colName) %>% toupper()
    plotSub <- paste0("n(genes) =", length(GENES), ", minTxLen =", minTxLength,
                      ", minReadsPerGene =", minReadsPerGene, ", smooth.sp =", smooth)
    Y_axis <- paste(FUN, colName)
    if(smooth){
        DF$smooth <- smooth.spline(DF$score)$y
        ggplot(DF, aes(x = bins, y = smooth)) + geom_line() + theme_classic() +
            ggtitle(plotTitle, plotSub) + ylab(Y_axis)
    }else{
        ggplot(DF, aes(x = bins, y = score)) + geom_line() + theme_classic() +
            ggtitle(plotTitle, plotSub) + ylab(Y_axis)
    }
}


tx_add_CtoTMR <- function(DT, minCov = 50, onNucs = c("C")){
    DT <- data.table::data.table(DT)
    tmp <- round(DT$`T` / (DT$`T` + DT$`C`), 6)
    tmp[DT$cov < minCov] <- NA
    tmp[!(DT$refSeq %in% onNucs)] <- NA
    tibble::add_column(DT, MR_CtoT = tmp)
}

# Adding rolling function to DT
tx_add_rollingMean <- function(DT, colName, winSize, newColName, fill = NA, align = "center", minCov = 21, nCores = 1){
    oNames <- colnames(DT)
    tmp <- mclapply(mc.cores = nCores, tx_split_DT(DT), function(x){
        RcppRoll::roll_mean(x[[colName]], n = winSize, fill = fill, align = align)
    }) %>% unlist() %>% unname()
    tmp[DT$cov < minCov] <- NA
    tibble::add_column(DT, tmp) %>% magrittr::set_names(c(oNames, newColName))
}

# shift column
tx_shift_geneWise <- function(DT, colToShift, direction, bp, nCores){
    if(!colToShift %in% colnames(DT)){
        stop("colToShift is not a column of DT")
    }
    DTL <- tx_split_DT(DT)
    if(direction == "upstream"){
        OUT <- mclapply(mc.cores = nCores, DTL, function(x){
            x[[colToShift]] <- c(tail(x[[colToShift]], -bp), rep(NA, bp))
            return(x)
        }) %>% tx_merge_DT()
    }else if(direction == "downstream"){
        OUT <- mclapply(mc.cores = nCores, DTL, function(x){
            x[[colToShift]] <- c(rep(NA, bp), head(x[[colToShift]], -bp))
            return(x)
        }) %>% tx_merge_DT()
    }else{
        stop("Input to argument 'direction' is not either 'downstream' or 'upstream'.")
    }
    return(OUT)
}

# ggseqlogo
tx_plot_ggseqlogo <- function(DT, annotCol, flankSize = 6, method = "bits", nCores = 1){
    ggOUT <- tx_get_flanksFromLogicAnnot(DT, annotCol = annotCol, valuesCol = "refSeq",
                                         flankSize = flankSize, nCores = nCores, minCov = 0) %>%
        apply(MARGIN =  1, function(x) paste(x, collapse = "")) %>%
        ggseqlogo::ggseqlogo(method = method)
    if(method == "bits"){
        ggOUT + ggplot2::ylim(0,2) + ggplot2::theme_minimal()
    }else{
        ggOUT + ggplot2::theme_minimal()
    }
}

tx_ggseqlogo <- tx_plot_ggseqlogo

# Combine lists of transcript reads processed by tx_reads
tx_combineTxReadsList <- function(txReadsList){
    tmp <- lapply(txReadsList, function(x) unlist(x)) %>%
        GenomicRanges::GRangesList() %>%
        unlist()
    split(tmp, seqnames(tmp))
}

# Linearized version
# Problem: will merge sequences between genes that are too close to each other.
tx_extractNeighSeq <- function(DT, var, upFlank, doFlank){
    Igenes <- unique(DT[DT[[var]],][["gene"]]) %>% as.character()
    DT <- DT[DT$gene %in% Igenes,]
    fullSeq <- paste0(DT$refSeq, collapse = "")
    str_sub(fullSeq, which(DT[[var]]) - upFlank, which(DT[[var]]) + doFlank)
}

# Function for cytidine persistence to Bisulphite treatment
tx_add_CtoTMR <- function(DT, minCov = 50, onNucs = c("C")){
    DT <- data.table::data.table(DT)
    tmp <- round(DT$`T` / (DT$`T` + DT$`C`), 6)
    tmp[DT$cov < minCov] <- NA
    tmp[!(DT$refSeq %in% onNucs)] <- NA
    tibble::add_column(DT, MR_CtoT = tmp)
}

# Function to extract the sequences of a transcriptome from a genome
tx_get_transcriptSeqs <- function(genome, geneAnnot, outFile = NULL, nCores = 1){
    check_GA_genome_chrCompat(geneAnnot, genome)
    CHRS <- unique(GenomicRanges::seqnames(geneAnnot))
    allSEQS <- parallel::mclapply(mc.cores = nCores, CHRS, function(iChr){
        subTXOME <- geneAnnot[which(as.logical(GenomicRanges::seqnames(geneAnnot) == iChr))]
        iBlocks <- IRanges::shift(S4Vectors::mcols(subTXOME)$blocks,
                                  IRanges::start(subTXOME) - 1)
        tmp <- stringr::str_sub(genome[[iChr]], start = IRanges::start(unlist(iBlocks)),
                                end = IRanges::end(unlist(iBlocks)))
        tmp3 <- lapply(split(tmp, rep(seq_along(iBlocks), times = sapply(iBlocks, length))),
                       function(x) paste(x, collapse = "")) %>%
            unlist() %>% Biostrings::DNAStringSet()
        tmp3[as.logical(GenomicRanges::strand(subTXOME) == "-")] <-
            Biostrings::reverseComplement(tmp3[as.logical(GenomicRanges::strand(subTXOME) == "-")])
        names(tmp3) <- subTXOME$name
        tmp3
    }) %>% do.call(what = "c")
    if(is.null(outFile)){
        allSEQS
    }else{
        Biostrings::writeXStringSet(allSEQS, filepath = outFile, format = "fasta")
    }
}

# Generate single-end FASTQ file
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
            DNAStringSet()
    }) %>% do.call(what = "c")
    names(metaTX$reads1) <- paste0("R1_", 1:length(metaTX$reads1))
    qualStr <- paste(rep("H", readLen), collapse = "")
    # Writing FASTQ
    Biostrings::writeXStringSet(x = metaTX$reads1,
                                quali = BStringSet(rep(qualStr, each = length(metaTX$reads1))),
                                filepath = fileName, format = "fastq", compress = TRUE)
}


# Homogenize genes in txDTs in a list
tx_unifyTxDTL <- function(txDTL, geneAnnot = NULL, genome = NULL, type = "intersect", nCores = 1){
    if(!all(Reduce(intersect, lapply(txDTL, colnames)) == Reduce(union, lapply(txDTL, colnames)))){
        stop("Column names of the elements of txDTL must be the same")
    }
    if(type == "intersection"){
        selGenes <- Reduce(intersect, lapply(txDTL, function(x) unique(x$gene)))
        new_txDTL <- mclapply(mc.cores = nCores, txDTL, function(x){
            x <- x[x$gene %in% selGenes]
            x[order(x$gene, x$txcoor)]
        })
        return(new_txDTL)
    }else if(type == "union"){
        if(is.null(geneAnnot)){stop("geneAnnot must be provided, as loaded with tx_load_bed()")}
        if(is.null(genome)){stop("geneAnnot must be provided, as loaded with tx_load_genome()")}
        selGenes <- Reduce(union, mclapply(mc.cores = nCores, txDTL, function(x) unique(x$gene)))
        selGA <- geneAnnot[geneAnnot$name %in% selGenes]
        new_txDTL <- mclapply(mc.cores = nCores, txDTL, function(x){
            x <- x[x$gene %in% selGenes]
            x <- x[order(x$gene, x$txcoor)]
            tx_complete_DT(DT = x, geneAnnot = selGA, nCores = nCores, genome = genome)
        })
    }else{stop("Argument 'type' must be either 'intersection' or 'union'")}
}

# shift column
tx_shift_geneWise <- function(DT, colToShift, direction, bp, nCores){
    if(!colToShift %in% colnames(DT)){
        stop("colToShift is not a column of DT")
    }
    DTL <- tx_split_DT(DT)
    if(direction == "upstream"){
        OUT <- mclapply(mc.cores = nCores, DTL, function(x){
            x[[colToShift]] <- c(tail(x[[colToShift]], -bp), rep(NA, bp))
            return(x)
        }) %>% tx_merge_DT()
    }else if(direction == "downstream"){
        OUT <- mclapply(mc.cores = nCores, DTL, function(x){
            x[[colToShift]] <- c(rep(NA, bp), head(x[[colToShift]], -bp))
            return(x)
        }) %>% tx_merge_DT()
    }else{
        stop("Input to argument 'direction' is not either 'downstream' or 'upstream'.")
    }
    return(OUT)
}

# TODO ####
# Make generic function to calculate the rate of nucleotide misincorporation
# set dividend and divisor nucleotides as well as a shift (1 downstream, 1 upstream, or 0)
# From brainSTORM

