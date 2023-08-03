# Plot meta gene by bins
tx_plot_metaGeneByBins <- function(DT, colName, nBins = 100, FUN = "mean", minTxLength = 300,
                                   minReadsPerGene = 100, smooth = TRUE, na.rm = TRUE){
    if(nBins >= minTxLength){stop("Number of bins most be smaller than minimum",
                                  "transcript length.")}
    if(!(colName %in% colnames(DT))){stop("colName is not a column in DT.")}
    nStarts <- DT[, sum(DT$start_5p), by = DT$gene][order(DT$V1, decreasing = T),]
    DT <- DT[DT$gene %in% nStarts$gene[nStarts$V1 >= minReadsPerGene]]
    geneLens <- tx_get_geneLengths(DT)
    DT <- DT[DT$gene %in% names(geneLens)[geneLens > minTxLength]]
    GENES <- as.character(unique(DT$gene))
    if(length(GENES) < 1){stop("No genes after filtering parameters")}
    meanCovBinned <- parallel::mclapply(seq(GENES), function(i){
        iGene <- GENES[i]
        tmpDT <- DT[DT$gene == iGene,]
        tmpDT$group <- tmpDT$txcoor %>% ggplot2::cut_number(n = nBins)
        out <- tapply(tmpDT[[colName]], tmpDT$group, FUN, na.rm = na.rm)
        out[is.nan(out)] <- NA
        out
    }) %>% do.call(what = "rbind") %>% apply(MARGIN = 2, FUN = FUN, na.rm = na.rm) %>%
        magrittr::set_names(NULL)
    DF <- data.frame(bins = seq(meanCovBinned) %>% as.numeric(),
                     score = meanCovBinned)
    plotTitle <- paste("METAGENE BY BINS -", FUN, colName) %>% toupper()
    plotSub <- paste0("n(genes) =", length(GENES), ", minTxLen =", minTxLength,
                      ", minReadsPerGene =", minReadsPerGene, ", smooth.sp =", smooth)
    Y_axis <- paste(FUN, colName)
    if(smooth){
        DF$smooth <- stats::smooth.spline(DF$score)$y
        ggplot2::ggplot(DF, ggplot2::aes(x = DF$bins, y = DF$smooth)) +
            ggplot2::geom_line() + ggplot2::theme_classic() +
            ggplot2::ggtitle(plotTitle, plotSub) + ggplot2::ylab(Y_axis)
    }else{
        ggplot2::ggplot(DF, ggplot2::aes(x = DF$bins, y = DF$score)) +
            ggplot2::geom_line() + ggplot2::theme_classic() +
            ggplot2::ggtitle(plotTitle, plotSub) + ggplot2::ylab(Y_axis)
    }
}

# tx_plot_metaGeneByBins <- function(DT, colName, nBins = 100, FUN = "mean", minTxLength = 300,
#                                    minReadsPerGene = 100, smooth = TRUE, na.rm = TRUE){
#     if(nBins >= minTxLength){stop("Number of bins most be smaller than minimum",
#                                   "transcript length.")}
#     if(!(colName %in% colnames(DT))){stop("colName is not a column in DT.")}
#     tapply(txDTn[["start_5p"]], txDT$gene, sum)
#     nStarts <- DT[, sum(start_5p), by = gene][order(V1, decreasing = T),]
#     DT <- DT[gene %in% nStarts$gene[nStarts$V1 >= minReadsPerGene]]
#     geneLens <- tx_get_geneLengths(DT)
#     DT <- DT[gene %in% names(geneLens)[geneLens > minTxLength]]
#     GENES <- as.character(unique(DT$gene))
#     if(length(GENES) < 1){stop("No genes after filtering parameters")}
#     meanCovBinned <- mclapply(seq(GENES), function(i){
#         iGene <- GENES[i]
#         tmpDT <- DT[gene == iGene,]
#         tmpDT$group <- tmpDT$txcoor %>% cut_number(n = nBins)
#         out <- tapply(tmpDT[[colName]], tmpDT$group, FUN, na.rm = na.rm)
#         out[is.nan(out)] <- NA
#         out
#     }) %>% do.call(what = "rbind") %>% apply(MARGIN = 2, FUN = FUN, na.rm = na.rm) %>% set_names(NULL)
#     DF <- data.frame(bins = seq(meanCovBinned) %>% as.numeric,
#                      score = meanCovBinned)
#     plotTitle <- paste("METAGENE BY BINS -", FUN, colName) %>% toupper()
#     plotSub <- paste0("n(genes) =", length(GENES), ", minTxLen =", minTxLength,
#                       ", minReadsPerGene =", minReadsPerGene, ", smooth.sp =", smooth)
#     Y_axis <- paste(FUN, colName)
#     if(smooth){
#         DF$smooth <- smooth.spline(DF$score)$y
#         ggplot(DF, aes(x = bins, y = smooth)) + geom_line() + theme_classic() +
#             ggtitle(plotTitle, plotSub) + ylab(Y_axis)
#     }else{
#         ggplot(DF, aes(x = bins, y = score)) + geom_line() + theme_classic() +
#             ggtitle(plotTitle, plotSub) + ylab(Y_axis)
#     }
# }

# load_all()

# Combine lists of transcript reads processed by tx_reads
tx_combineTxReadsList <- function(txReadsList){
    tmp <- lapply(txReadsList, function(x) unlist(x)) %>%
        GenomicRanges::GRangesList() %>%
        unlist()
    split(tmp, GenomicAlignments::seqnames(tmp))
}

# Ideas for functions ##########################################################

# # Apply a function in a binned manner to a data.table column
# tx_DT_bin_function <- function(DT, col, bins = 100, fun = mean){
#     tapply(X = DT[[col]], INDEX = cut_interval(1:nrow(DT), bins), FUN = fun) %>% set_names(c())
# }

# # Apply a function in a binned manner to a vector
# tx_bin_function <- function(x, bins = 100, fun = mean){
#     tapply(X = x, INDEX = cut_interval(1:length(x), bins), FUN = fun) %>% set_names(c())
# }


#' Generate paired-end FASTQ file
#'
#' Simulates a paired-end FASTQ files from a genome and a gene annotation.
#' The distribution of reads is randomly selected following a negative binomial
#' distribution.
#'
#' @param genome list. The full reference genome sequences, as loaded by
#' \code{\link{tx_load_genome}}() or prepackaged by BSgenome, see ?BSgenome::available.genomes
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().
#' @param readLen_1 integer. Length of simulated reads1
#' @param readLen_2 integer. Length of simulated reads2
#' @param insertSize  integer. Length of simulated insert (gap between reads)
#' @param libSize integer. Size of simulated FASTQ file
#' @param TPM numeric. Proportion of reads per gene in gene annotation.
#' Overrides negative binomial distribution generation.
#' @param NB_r numeric. Target r dispersion parameter. Bigger values of r tend
#' to generate a normal distribution. See ?stats::rnbinom()
#' @param NB_mu numeric. Target mean of the resulting distribution
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability not available in Windows OS.
#' @param fileName_r1  character. Name of output file reads1.
#' @param fileName_r2  character. Name of output file reads2.
#'
#' @return Writes FASTQ files fileName_r1 and fileName_r2
#'
tx_generatePairedEndFASTQ <- function(genome, geneAnnot, readLen_1, readLen_2,
                                      insertSize, libSize, fileName_r1,
                                      fileName_r2, TPM = NULL, NB_r = 5, NB_mu = 500, nCores){
    # filter transcripts by size
    totReadLen <- (readLen_1 + readLen_2 + insertSize)
    txOME_seqs <- tx_get_transcriptSeqs(genome = genome, geneAnnot = geneAnnot,
                                        nCores = nCores)
    txOME_seqs <- txOME_seqs[BiocGenerics::width(txOME_seqs) >= totReadLen]
    if(length(txOME_seqs) != length(geneAnnot)){
        stop("Some genes in gene annotation are shorter than total read length\n")
    }
    if(!is.null(TPM)){
        if(length(TPM) != length(geneAnnot)){
            stop("If TPM is provided, its length must be equal to geneAnnot's length")
        }
    }
    # Random starts
    metaTX <- list(seqs = txOME_seqs, width = BiocGenerics::width(txOME_seqs))
    if(is.null(TPM)){
        metaTX$TPM <- stats::rnbinom(n = length(txOME_seqs), size = NB_r, mu = NB_mu)
        metaTX$TPM <- round(metaTX$TPM * (libSize / sum(metaTX$TPM)))
    }else{
        metaTX$TPM <- TPM
        metaTX$TPM <- round(metaTX$TPM * (libSize / sum(metaTX$TPM)))
    }
    metaTX$range1 <- metaTX$width - (totReadLen) + 1
    metaTX$starts_1 <- lapply(seq_along(metaTX$TPM), function(i){
        sample(seq(1, metaTX$range1[i]), size = metaTX$TPM[i], replace = TRUE)
    })
    #Extract sequences
    metaTX$reads1 <- lapply(seq_along(metaTX$seqs), function(i){
        stringr::str_sub(metaTX$seqs[i], metaTX$starts_1[[i]],
                         metaTX$starts_1[[i]] + readLen_1 - 1) %>%
            Biostrings::DNAStringSet()
    }) %>% do.call(what = "c")
    metaTX$reads2 <- lapply(seq_along(metaTX$seqs), function(i){
        stringr::str_sub(metaTX$seqs[i], metaTX$starts_1[[i]] + readLen_1 + insertSize,
                         metaTX$starts_1[[i]] + totReadLen - 1) %>%
            Biostrings::DNAStringSet() %>% Biostrings::reverseComplement()
    }) %>% do.call(what = "c")
    names(metaTX$reads1) <- paste0("R1_", 1:length(metaTX$reads1))
    names(metaTX$reads2) <- paste0("R2_", 1:length(metaTX$reads2))
    qualStr_1 <- paste(rep("H", readLen_1), collapse = "")
    qualStr_2 <- paste(rep("H", readLen_2), collapse = "")
    # Writing FASTQS
    Biostrings::writeXStringSet(x = metaTX$reads1,
                                quali = Biostrings::BStringSet(
                                    rep(qualStr_1, each = length(metaTX$reads1))),
                                filepath = fileName_r1,
                                format = "fastq",
                                compress = TRUE)
    Biostrings::writeXStringSet(x = metaTX$reads2,
                                quali = Biostrings::BStringSet(
                                    rep(qualStr_2, each = length(metaTX$reads2))),
                                filepath = fileName_r2,
                                format = "fastq",
                                compress = TRUE)
}

#' Generate single-end FASTQ file
#'
#' Simulates a single-end FASTQ file from a genome and a gene annotation.
#' The distribution of reads is randomly selected following a negative binomial
#' distribution.
#'
#' @param genome list. The full reference genome sequences, as loaded by
#' \code{\link{tx_load_genome}}() or prepackaged by BSgenome, see ?BSgenome::available.genomes
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().
#' @param readLen integer. Length of simulated reads
#' @param libSize integer. Size of simulated FASTQ file
#' @param fileName character. Name of output file.
#' @param TPM numeric. Proportion of reads per gene in gene annotation
#' @param NB_r numeric. Target r dispersion parameter. Bigger values of r tend
#' to generate a normal distribution. See ?stats::rnbinom()
#' @param NB_mu numeric. Target mean of the resulting distribution
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability not available in Windows OS.
#'
#' @return Writes FASTQ fileName
tx_generateSingleEndFASTQ <- function(genome, geneAnnot, readLen, libSize,
                                      fileName, TPM = NULL, NB_r = 5, NB_mu = 500, nCores){
    # filter transcripts by size
    txOME_seqs <- tx_get_transcriptSeqs(genome = genome, geneAnnot = geneAnnot,
                                        nCores = nCores)
    txOME_seqs <- txOME_seqs[BiocGenerics::width(txOME_seqs) >= readLen]
    if(length(txOME_seqs) != length(geneAnnot)){
        stop("Some genes in gene annotation are shorter than total read length\n")
    }
    if(!is.null(TPM)){
        if(length(TPM) != length(geneAnnot)){
            stop("If TPM is provided, its length must be equal to geneAnnot's length")
        }
    }
    # Random starts
    metaTX <- list(seqs = txOME_seqs, width = BiocGenerics::width(txOME_seqs))
    if(is.null(TPM)){
        metaTX$TPM <- stats::rnbinom(n = length(txOME_seqs), size = NB_r, mu = NB_mu)
        metaTX$TPM <- round(metaTX$TPM * (libSize / sum(metaTX$TPM)))
    }else{
        metaTX$TPM <- TPM
        metaTX$TPM <- round(metaTX$TPM * (libSize / sum(metaTX$TPM)))
    }

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
                                quali = Biostrings::BStringSet(
                                    rep(qualStr, each = length(metaTX$reads1))),
                                filepath = fileName, format = "fastq",
                                compress = TRUE)
}
