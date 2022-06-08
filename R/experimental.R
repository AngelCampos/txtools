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
    }) %>% do.call(what = rbind) %>% apply(MARGIN = 2, FUN = FUN, na.rm = na.rm) %>%
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
        ggplot2::ggplot(DF, ggplot2::aes(x = DF$bins, y = DF$score)) + ggplot2::geom_line() + ggplot2::theme_classic() +
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
#     }) %>% do.call(what = rbind) %>% apply(MARGIN = 2, FUN = FUN, na.rm = na.rm) %>% set_names(NULL)
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

# Plot one numeric variable
tx_plot_oneNumeric <- function(DTL, gene, txRange = 1:nrow(DTL[[1]]), columnName,
                               show_yLabels = TRUE, bar_border = TRUE,
                               showLegend = TRUE, col = "#c2c2c2"){
    if(all(methods::is(DTL) %in% methods::is(list()))){
        lapply(DTL, function(x) check_refSeq(x))
        DTL <- lapply(DTL, function(x) check_DT(x))
        if(!all(unlist(lapply(DTL, function(x) gene %in% x$gene)))){stop("gene not found in DT object")}
        if(is.null(names(DTL))){names(DTL) <- seq_along(DTL)}
        tmpData <- lapply(seq_along(DTL), function(i){
            DT <- DTL[[i]]
            DT <- DT[DT$gene == gene, ]
            DT <- DT[DT$txcoor %in% txRange,]
            DT$pos <- paste(DT$txcoor, DT$refSeq, sep = "-")
            DT$pos <- factor(DT$pos, levels = DT$pos)
            tmpData <- tidyr::pivot_longer(DT, cols = tidyr::all_of(columnName),
                                           values_to = "counts",
                                           names_to = "coverage") %>%
                data.table::data.table()
            tmpData$coverage <- factor(tmpData$coverage, levels = columnName)
            tmpData$sample <- names(DTL)[i]
            tmpData
        }) %>% do.call(what = rbind)
    }else{
        DT <- DTL
        DT <- check_refSeq(DT)
        DT <- check_DT(DT)
        DT <- as.data.frame(DT)
        if(!(gene %in% DT$gene)){stop("gene not found in DT object")}
        DT <- DT[DT$gene == gene, ]
        DT <- DT[DT$txcoor %in% txRange,]
        DT$pos <- paste(DT$txcoor, DT$refSeq, sep = "-")
        DT$pos <- factor(DT$pos, levels = DT$pos)
        tmpData <- tidyr::pivot_longer(DT, cols = tidyr::all_of(columnName),
                                       values_to = "counts",
                                       names_to = "coverage") %>%
            data.table::data.table()
        tmpData$coverage <- factor(tmpData$coverage, levels = columnName)
    }
    tmpGG <- ggplot2::ggplot(tmpData,
                             ggplot2::aes(x = tmpData$pos,
                                          y = tmpData$counts,
                                          fill = tmpData$coverage)) +
        ggplot2::theme_minimal() +
        ggplot2::scale_fill_manual(values = col) +
        ggplot2::ylab(columnName) + ggplot2::xlab("Transcriptome coordinate") +
        ggplot2::theme(legend.position="none") +
        ggplot2::guides(fill = ggplot2::guide_legend(nrow=1, byrow=TRUE, title = "")) +
        ggplot2::ggtitle(gene)
    if(bar_border){
        tmpGG <- tmpGG + ggplot2::geom_bar(stat = "identity",
                                           colour = "black",
                                           size = 0.3)
    }else{
        tmpGG <- tmpGG + ggplot2::geom_bar(stat = "identity")
    }
    if(show_yLabels){
        if(all(methods::is(DTL) == methods::is(list()))){
            tmpDT <- DTL[[1]]
            nucCols <- txBrowser_pal()(6)[-1:-2][as.numeric(
                factor(tmpDT[tmpDT$gene == gene & tmpDT$txcoor %in% txRange, ]$refSeq))]
        }else{
            nucCols <- txBrowser_pal()(6)[-1:-2][as.numeric(factor(DTL$refSeq))]
        }
        tmpGG <- suppressWarnings(
            tmpGG +
                ggplot2::theme(axis.text.x =
                                   ggplot2::element_text(angle = 90,
                                                         hjust = 1,
                                                         vjust = 0.5,
                                                         colour = nucCols,
                                                         face = "bold")))
    }else{
        tmpGG <- tmpGG + ggplot2::theme(axis.text.x = ggplot2::element_blank())
    }
    if(!showLegend){
        tmpGG <- tmpGG + ggplot2::theme(legend.position="none")
    }
    if(all(methods::is(DTL) == methods::is(list()))){
        tmpGG <- tmpGG + ggplot2::facet_wrap(sample~.)
    }else{
        tmpGG
    }
}

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
