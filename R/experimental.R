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

# Combine lists of transcript reads processed by tx_reads
tx_combineTxReadsList <- function(txReadsList){
    tmp <- lapply(txReadsList, function(x) unlist(x)) %>%
        GenomicRanges::GRangesList() %>%
        unlist()
    split(tmp, seqnames(tmp))
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

# Make generic function to calculate the rate of nucleotide misincorporation
# set dividend and divisor nucleotides as well as a shift (1 downstream, 1 upstream, or 0)
# From brainSTORM


check_sameGenesInDTL <- function(DTL){
    tmpL <- lapply(DTL, function(x) base::unique(x$gene))
    tmpA <- base::Reduce(x = tmpL, union)
    tmpB <- base::Reduce(x = tmpL, intersect)
    all(tmpA %in% tmpB)
}

check_DThasCol <- function(DT, colName){
    DT <- check_DT(DT)
    if(!colName %in% colnames(DT)){stop("DT does not contain '", colName, "' column.")}
}

# Meta gene #####

annot_CDSsta_DTL <- Vectorize(FUN = function(DT, CDS_start){
        DT$CDS_start[DT$gencoor == CDS_start[CDS_start$gene == DT$gene[1],]$gencoor] <- TRUE
        DT}, vectorize.args = "DT", SIMPLIFY = FALSE)

annot_CDSend_DTL <- Vectorize(FUN = function(DT, CDS_end){
    DT$CDS_end[DT$gencoor == CDS_end[CDS_end$gene == DT$gene[1],]$gencoor] <- TRUE
    DT}, vectorize.args = "DT", SIMPLIFY = FALSE)

check_GA_txDT_compat <- function(DT, geneAnnot){
    check_DThasCol(DT, "gene")
    if(!all(as.character(unique(DT$gene)) %in% geneAnnot$name)){
        stop("Not all genes in DT are contained in geneAnnot")
    }
}

# Potential functions ##########################################################

# # Apply a function in a binned manner to a data.table column
# tx_DT_bin_function <- function(DT, col, bins = 100, fun = mean){
#     tapply(X = DT[[col]], INDEX = cut_interval(1:nrow(DT), bins), FUN = fun) %>% set_names(c())
# }
#
# # Apply a function in a binned manner to a vector
# tx_bin_function <- function(x, bins = 100, fun = mean){
#     tapply(X = x, INDEX = cut_interval(1:length(x), bins), FUN = fun) %>% set_names(c())
# }
