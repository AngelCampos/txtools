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

# # Adding rolling function to DT
tx_add_rollingMean <- function(DT, colName, winSize, newColName = NULL, fill = NA, align = "center", minCov = 21, nCores = 1){
    if(is.null(newColName)){newColName <- paste(colName , "rollMean", sep = "_" )}
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

# Combine lists of transcript reads processed by tx_reads
tx_combineTxReadsList <- function(txReadsList){
    tmp <- lapply(txReadsList, function(x) unlist(x)) %>%
        GenomicRanges::GRangesList() %>%
        unlist()
    split(tmp, seqnames(tmp))
}

# Linearized version tx_get_neigh...
tx_get_flankSequence <- function(DT, logi_col, upFlank, doFlank, addNames = FALSE){
    if(!is.logical(DT[[logi_col]])){stop("column ", logi_col, " must be of class 'logical'.")}
    if(!"refSeq" %in% colnames(DT)){stop("refSeq must be in DT. You can add it with tx_add_refSeq().")}
    if(!"gene" %in% colnames(DT)){stop("gene column must be in DT.")}
    Igenes <- as.character(unique(DT[DT[[logi_col]],][["gene"]]))
    DT <- DT[DT$gene %in% Igenes,]
    tmpSeq <- split(DT$refSeq, forcats::fct_drop(DT$gene))
    tmpVar <- split(DT[[logi_col]], forcats::fct_drop(DT$gene))
    spacer <- rep(".", max(upFlank, doFlank))
    spacerVar <- rep(NA, max(upFlank, doFlank))
    fullSeq <- c(spacer, lapply(tmpSeq, function(x) (c(x, spacer))) %>%
                     do.call(what = "c")) %>% paste(collapse = "")
    fullVar <- c(spacerVar, lapply(tmpVar, function(x) (c(x, spacerVar))) %>% do.call(what = "c"))
    tmpSeq <- stringr::str_sub(fullSeq, which(fullVar) - upFlank, which(fullVar) + doFlank)
    if(any(grepl(tmpSeq, pattern = "\\."))){
        warning("Some sequences reached the end of transcript, a '.'",
                "was added in place, which may affect downstream results.")
    }
    if(addNames){
        names(tmpSeq) <- paste(DT$gene[DT[[logi_col]]], DT$txcoor[DT[[logi_col]]], sep = ":")
    }
    return(tmpSeq)
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
tx_unifyTxDTL <- function(txDTL, geneAnnot = NULL, genome = NULL, type = "intersection", nCores = 1){
    if(!all(Reduce(intersect, lapply(txDTL, colnames)) == Reduce(union, lapply(txDTL, colnames)))){
        stop("Column names of the elements of txDTL must be the same")
    }
    if(type == "intersection"){
        selGenes <- Reduce(intersect, lapply(txDTL, function(x) unique(x$gene)))
        new_txDTL <- parallel::mclapply(mc.cores = nCores, txDTL, function(x){
            x <- x[x$gene %in% selGenes,]
            x[order(x$gene, x$txcoor),]
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
tx_shift_geneWise <- function(DT, colToShift, direction, bp, nCores = 1){
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

tx_test_ttest <- function(DTL, cont_var, test_groups, test_na.rm = FALSE, ...){
    sapply(DTL, function(x) check_DThasCol(x, cont_var))
    #check if needed to unify
    if(!check_sameGenesInDTL(DTL)){
        DTL <- tx_unifyTxDTL(DTL, type = "intersection")
    }
    test_mat <- do.call(what = "cbind", lapply(DTL, function(x) x[[cont_var]]))
    rowTtests <- genefilter::rowttests(x = test_mat, fac = test_groups, na.rm = test_na.rm, ...)
    if("refSeq" %in% colnames(DTL[[1]])){
        txRES <- DTL[[1]][, c("chr", "gencoor", "strand", "gene", "txcoor", "refSeq")]
    }else{
        txRES <- DTL[[1]][,c("chr", "gencoor", "strand", "gene", "txcoor")]
    }
    txRES <- cbind(txRES, rowTtests)
}

annot_CDSsta_DTL <- Vectorize(FUN = function(DT){
        DT$CDS_start[DT$gencoor == CDS_start[CDS_start$gene == DT$gene[1],]$gencoor] <- TRUE
        DT}, vectorize.args = "DT", SIMPLIFY = FALSE)

annot_CDSend_DTL <- Vectorize(FUN = function(DT){
    DT$CDS_end[DT$gencoor == CDS_end[CDS_end$gene == DT$gene[1],]$gencoor] <- TRUE
    DT}, vectorize.args = "DT", SIMPLIFY = FALSE)

check_GA_txDT_compat <- function(DT, geneAnnot){
    check_DThasCol(DT, "gene")
    if(!all(as.character(unique(DT$gene)) %in% geneAnnot$name)){
        stop("Not all genes in DT are contained in geneAnnot")
    }
}

tx_get_metageneAtCDS <- function(txDT, geneAnnotation, colVars, CDS_align, upFlank, doFlank){
    check_GA_txDT_compat(txDT, geneAnnotation)
    invisible(sapply(c("gencoor", "strand", "gene"), function(x) check_DThasCol(txDT, x)))
    geneAnnotation <- geneAnnotation[geneAnnotation$name %in% unique(txDT$gene)]
    NCG <- geneAnnotation[GenomicRanges::width(geneAnnotation$thick) == 0]
    CG <- geneAnnotation[GenomicRanges::width(geneAnnotation$thick) > 0]
    if(length(NCG)>0){warning(length(NCG), " non-coding genes where ommitted from analysis.")}
    if(sum(GenomicRanges::strand(CG) == "*") > 0){
        stop("Genes with no set strand (*) are not allowed in geneAnnot.")
    }
    pos_CG <- as.factor(GenomicRanges::strand(CG)) == "+"
    neg_CG <- as.factor(GenomicRanges::strand(CG)) == "-"
    if(CDS_align == "start"){
        CDS_start <- rbind(data.frame(gene = CG[pos_CG]$name,
                                      gencoor = IRanges::start(CG[pos_CG]$thick)),
                           data.frame(gene = CG[neg_CG]$name,
                                      gencoor = IRanges::end(CG[neg_CG]$thick)))
        txDT$CDS_start <- FALSE
        txDT <- tx_split_DT(txDT) %>% annot_CDSsta_DTL() %>% tx_merge_DT()
        tmpFlanks <- lapply(colVars, function(colVar){
            tx_get_flanksFromLogicAnnot(DT = txDT,
                                        logi_col = "CDS_start",
                                        upFlank = upFlank,
                                        doFlank = doFlank,
                                        values_col = colVar,
                                        addRowNames = TRUE)})
    }else if(CDS_align == "end"){
        CDS_end <- rbind(data.frame(gene = CG[pos_CG]$name,
                                    gencoor = IRanges::end(CG[pos_CG]$thick)),
                         data.frame(gene = CG[neg_CG]$name,
                                    gencoor = IRanges::start(CG[neg_CG]$thick)))
        
        txDT$CDS_end <- FALSE
        txDT <- tx_split_DT(txDT) %>% annot_CDSend_DTL() %>% tx_merge_DT()
        tmpFlanks <- lapply(colVars, function(colVar){
            tx_get_flanksFromLogicAnnot(DT = txDT,
                                        logi_col = "CDS_end",
                                        upFlank = upFlank,
                                        doFlank = doFlank,
                                        values_col = colVar,
                                        addRowNames = TRUE)
        }) 
    }else{stop("CDS_align should be either 'start' or 'end'.")}
    return(tmpFlanks %>% magrittr::set_names(colVars))
}

tx_plot_metageneAtCDS <- function(txDT, geneAnnotation, colVars, CDS_align, upFlank,
                                  doFlank, summ_fun = "sum", roll_fun = "sum", roll_n = 100,
                                  roll_align = "center", roll_fill = NA, smooth = TRUE, spar  = 0.3,
                                  na.rm = TRUE, normalize = TRUE, tick_by = NULL){
    tmpO <- tx_get_metageneAtCDS(txDT = txDT, geneAnnotation = geneAnnotation, colVars = colVars,
                                 CDS_align = CDS_align, upFlank = upFlank, doFlank = doFlank)
    tmpDF <- lapply(names(tmpO), function(x){
        if(summ_fun == "sum"){
            tmp2 <- colSums(tmpO[[x]], na.rm = na.rm)
        }else if(summ_fun == "mean"){
            tmp2 <- colMeans(tmpO[[x]], na.rm = na.rm)
        }else{stop("Argument 'summ_fun' has to be either sum or mean")}
        if(is.null(roll_fun)){
            tmp3 <- tmp2
        }else if(roll_fun == "sum"){
            tmp3 <- RcppRoll::roll_sum(tmp2, n = roll_n, align = roll_align, fill = roll_fill, na.rm = na.rm)
        }else if(roll_fun == "mean"){
            tmp3 <- RcppRoll::roll_mean(tmp2, n = roll_n, align = roll_align, fill = roll_fill, na.rm = na.rm)
        }else{stop("Argument 'roll_fun' has to be either sum or mean")}
        tmpDF <- data.frame(value = tmp3, position = factor(names(tmp2), levels = names(tmp2)), group = x)
        if(smooth){
            tmpDF[!is.na(tmpDF$value), ]$value <- stats::smooth.spline(tmpDF[!is.na(tmpDF$value), ]$value, spar = spar)$y
        }
        if(normalize){
            tmpDF$value <- (tmpDF$value / sum(tmpDF$value, na.rm = TRUE)) * 100
        }
        tmpDF
    }) %>% do.call(what = "rbind") %>% data.table::data.table()
    if(is.null(tick_by)){
        tick_by <- upFlank / 2
    }
    tmpGG <- ggplot(tmpDF, aes(x = position, y = value, group = group, colour = group)) +
        geom_line() + 
        scale_x_discrete(limits = unique(tmpDF$position),
                         breaks = unique(tmpDF$position)[seq(1, length(unique(tmpDF$position)), by = tick_by)]) +
        theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
    if(CDS_align == "end"){
        tmpGG <- tmpGG + ggplot2::geom_vline(xintercept = "CDS_end" , col = "black", linetype="dashed") + 
            ggtitle("Metagene aligned at CDS_end")
    }else if(CDS_align == "start"){
        tmpGG <- tmpGG + ggplot2::geom_vline(xintercept = "CDS_start" , col = "black", linetype="dashed") + 
            ggtitle("Metagene aligned at CDS_start")
    }
    tmpGG
}

