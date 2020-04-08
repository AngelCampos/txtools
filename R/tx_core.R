#' Read paired end bam file by yield size
#'
#' Reads a file in BAM format by blocks of lines equal to a yield size, either
#' automatically calculated or specified by the user, and loads it as a
#' GenomicAlignments object.
#'
#' @param file character. Path to the file to read
#' @param yieldSize numeric. Number of reads to be processed at a time
#' @param scanFlag integer. Flag used to filter reads. Check ?scanBamFlag()
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
#' @author M.A. Garcia-Campos
#'
#' @examples
#' # Loading  in-package BAM file
#' bamFile <- system.file("extdata", "example_hg19.bam", package = "txtools")
#' tx_load_bam(bamFile, loadSeq = TRUE, verbose = TRUE)
tx_load_bam <- function(file,
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
    if(verbose){
        cat("Loading BAM file \n")
        pb <- utils::txtProgressBar(style = 3)
    }
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
                                                   param = Rsamtools::ScanBamParam(flag = scanFlag,
                                                                                   what = "seq"))
        }else if(loadSeq == F){
            GenomicAlignments::readGAlignmentPairs(BAMFILE, use.names = TRUE,
                                                   param = Rsamtools::ScanBamParam(flag = scanFlag))
        }else(stop("loadSeq parameter must be logical either FALSE or TRUE"))
    }) %>% do.call(what = c)
    close(BAMFILE)
    if(verbose){close(pb)}
    bamData <- list(GAligns = reads,
                    dumpedAmbigPairs = GenomicAlignments::getDumpedAlignments())
    if(verbose){cat(" \n")}
    if(verbose){
        cat(length(bamData$GAligns), "reads succesfully loaded \n")
        cat("Dumped reads due to ambiguous pairs:", length(bamData$dumpedAmbigPairs), "\n")
    }

    if(recoverDumpedAligns == F){
        return(reads)
    }else{
        return(bamData)
    }
}

#' Load gene models from bed-12 and bed-6 files
#'
#' Reads and loads bed files checking for compliance with bed file structure
#'
#' @param bedfile character. Gene annotation file in bed-6 or bed-12 format
#'
#' @return
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_load_bed <- function(bedfile){
    tmp <- plyranges::read_bed(bedfile)
    if(length(S4Vectors::mcols(tmp)) == 2){
        tmp$itemRGgb <- NA
        tmp$thick <- tmp@ranges
        tmp$blocks <- IRanges::IRanges(start = 1,
                                       end = IRanges::end(tmp) -
                                           IRanges::start(tmp) + 1)
    }
    return(tmp)
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
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_reads <- function(reads, bedR, overlapType = "within", minReads = 50,
                     withSeq = F, verbose = T){
    if(class(reads) != "GAlignmentPairs"){
        stop("reads argument should be of class GAlignmentPairs \n")
    }
    if(class(bedR) != "GRanges"){
        stop("bedR argument should be of class GRanges \n")
    }
    if(verbose){
        cat("Processing", length(reads), "reads, using", length(bedR),
            "gene models \n")
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
    OUT <- OUT[lapply(OUT, length) %>% unlist %>%
                   magrittr::is_greater_than(minReads)] %>%
        GenomicRanges::GRangesList()
    if(verbose){
        cat("Output contains:", lapply(OUT, names) %>% unlist %>% unique %>% length,
            "unique reads in", length(OUT), "gene models \n")
    }
    return(OUT)
}

#' Merge paired-end reads and assign them to gene models (Multicore)
#'
#' @param reads
#' @param bedR
#' @param nCores
#' @param overlapType
#' @param minReads
#' @param withSeq
#' @param verbose
#'
#' @return
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_reads_mc <- function(reads, bedR, nCores, overlapType = "within",
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
        cat("Processing", length(reads), "paired-end reads, using", length(bedR),
            "gene models \n")
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
    OUT <- OUT[lapply(OUT, length) %>% unlist %>% magrittr::is_greater_than(minReads)] %>%
        GenomicRanges::GRangesList()
    if(verbose){
        cat("Output contains:", lapply(OUT, names) %>% unlist %>% unique %>% length,
            "unique reads in", length(OUT), "gene models \n")
    }
    return(OUT)
}

#' Filter ranges by a maximum width
#'
#' @param x GenomicAlignment
#' @param thr
#'
#' @return
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_filter_max_width <- function(x, thr){
    tmp <- GenomicAlignments::width(x) %>% magrittr::is_weakly_less_than(thr)
    lapply(seq(1, length(x)), function(i){
        if(all(tmp[[i]])){
            x[[i]]
        }else{
            x[[i]][tmp[[i]]]
        }
    }) %>% GenomicRanges::GRangesList() %>% magrittr::set_names(names(x))
}

#' Calculate coverage table cov, ends, starts for all gene models
#'
#' @param x
#'
#' @return
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_coverageTab <- function(x){
    if(class(x) != "CompressedGRangesList"){
        stop("x must be of class SimpleGRangesList")
    }
    if(names(x) %>% duplicated() %>% sum() %>% magrittr::is_greater_than(0)){
        stop("List contains duplicated gene models")
    }
    cov <- tx_coverage(x) %>% lapply(as.vector)
    sta <- GenomicRanges::start(x) %>% lapply(table) %>% magrittr::set_names(names(x))
    end <- GenomicRanges::end(x) %>% lapply(table) %>% magrittr::set_names(names(x))
    len <- GenomeInfoDb::seqlengths(x)
    OUT <- lapply(1:length(x), function(i){
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

#' Calculate coverage table (cov, ends, starts) for all gene models (Multicore)
#'
#' @param x
#' @param nCores
#'
#' @return
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_coverageTab_mc <- function(x, nCores){
    if(class(x) != "CompressedGRangesList"){
        stop("x must be of class SimpleGRangesList")
    }
    if(names(x) %>% duplicated() %>% sum() %>% magrittr::is_greater_than(0)){
        stop("List contains duplicated gene models")
    }
    cov <- tx_coverage(x) %>% mclapply(as.vector, mc.cores = nCores)
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
#' @param x
#' @param simplify_IUPAC
#'
#' @return
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_nucFreqTab <- function(x, simplify_IUPAC = "not"){
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

#' Table with genomic and transcriptomic coordinates
#'
#' @param x
#' @param geneAnnot_GR
#'
#' @return
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_genCoorTab <- function(x, geneAnnot_GR){
    if(all(names(x) %in% geneAnnot_GR$name)){
        lapply(names(x), function(iGene){
            tmp2 <- geneAnnot_GR[which(geneAnnot_GR$name == iGene)]
            tmp3 <- c(GenomicAlignments::seqnames(tmp2),
                      GenomicRanges::strand(tmp2)) %>% as.character() %>% c(iGene)
            tmpDT <- rep(tmp3, GenomeInfoDb::seqlengths(x[[iGene]])[iGene]) %>%
                matrix(ncol = 3, byrow = T) %>%
                cbind(exonBlockGen(iGene, geneAnnot_GR)) %>%
                cbind(seq(1, GenomeInfoDb::seqlengths(x[[iGene]])[iGene])) %>%
                .[,c(1,4,2,3,5)] %>% data.table::data.table() %>%
                magrittr::set_colnames(c("chr", "gencoor", "strand", "gene", "txcoor"))
            tmpDT$chr <- as.factor(tmpDT$chr)
            tmpDT$gencoor <- as.integer(tmpDT$gencoor)
            tmpDT$strand <- as.factor(tmpDT$strand)
            tmpDT$gene <- as.factor(tmpDT$gene)
            tmpDT$txcoor <- as.integer(tmpDT$txcoor)
            return(tmpDT)
        }) %>% magrittr::set_names(names(x))
    }else{
        stop("Names of x are not contained in geneAnnot_GR$name")
    }
}

#' Summarized Coverage data.table
#'
#' @param x
#' @param geneAnnot
#'
#' @return data.table
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_coverageDT <- function(x, geneAnnot){
    hlp_cbind2Tabs(tx_genCoorTab(x, geneAnnot),
                   tx_coverageTab(x))
}

#' Summarized Nucleotide Frequency data.table
#'
#' @param x GenomicRanges. Trancriptomic converted RNA-seq reads
#' @param geneAnnot
#' @param simplify_IUPAC string. Available options are :
#' \itemize{
#' \item "not": Will output the complete nucleotide frequency table including
#' ambiguous reads using the IUPAC ambiguity code. Check: ?Biostrings::IUPAC_CODE_MAP
#' \item "splitForceInt": Will force an integers split in which ambiguous codes
#' will be split and assigned half the frequency into their respective nucleotides,
#' if the frequency is an odd number the uneven count will be assigned as "N".
#' \item "splitHalf": Ambiguous nucleotide frequencies will be split in half to
#' their corresponding nucleotides, in cases where frequency is odd creating
#' non-integer frequencies.
#' }
#'
#' @return data.table
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_nucFreqDT <- function(x, geneAnnot, simplify_IUPAC = "splitForceInt"){
    hlp_cbind2Tabs(tx_genCoorTab(x, geneAnnot),
                   tx_nucFreqTab(x, simplify_IUPAC))
}

#' Summarized Coverage & Nucleotide Frequency data.table
#'
#' @param x
#' @param geneAnnot
#' @param simplify_IUPAC string. Available options are :
#' \itemize{
#' \item "not": Will output the complete nucleotide frequency table including
#' ambiguous reads using the IUPAC ambiguity code. Check: ?Biostrings::IUPAC_CODE_MAP
#' \item "splitForceInt": Will force an integers split in which ambiguous codes
#' will be split and assigned half the frequency into their respective nucleotides,
#' if the frequency is an odd number the uneven count will be assigned as "N".
#' \item "splitHalf": Ambiguous nucleotide frequencies will be split in half to
#' their corresponding nucleotides, in cases where frequency is odd creating
#' non-integer frequencies.
#' }
#'
#' @return data.table
#' @export
#'
#' @examples
tx_covNucFreqDT <- function(x, geneAnnot, simplify_IUPAC = "splitForceInt"){
    hlp_cbind3Tabs(tx_genCoorTab(txReads, geneAnnot),
                   tx_coverageTab(txReads),
                   tx_nucFreqTab(txReads, simplify_IUPAC))
}

#' Quantifies reads by gene model from a tx_reads list
#'
#' @param x GenomicRanges.
#'
#' @return
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_counts <- function(x){
    suppressWarnings(unlist(x)) %>% GenomeInfoDb::seqnames() %>% table()
}

#' Add reference sequence
#'
#' @param DT data.table. A summarized data.table object. See tx_coverageDT(),
#' tx_nucFreq() and tx_covNucFreqDT() functions.
#' @param fastaGenome list. The full reference genome sequences, as prepackaged
#' by BSgenome. See ?BSgenome::available.genomes()
#' @param geneAnnot GRanges. Gene annotation loaded as a GenomicRanges object,
#'  see tx_load_bed().
#'
#' @return
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_addRefSeqDT <- function(DT, fastaGenome, geneAnnot){
    if(class(DT) != "list"){
        stop("DT must be of class list")
    }
    allDT <- lapply(DT, function(x){
        "data.table" %in% class(x)
    }) %>% unlist %>% mean %>% magrittr::is_less_than(1)
    if(allDT){
        stop("All elements of DT must be of data.table class")
    }
    lapply(seq(DT), function(i){
        iGene <- as.character(unique(DT[[i]]$gene))
        if(length(iGene) != 1){stop("Each element in DT must be a data.table
                                    pertaining only one gene")}
        iChr <- as.character(DT[[i]]$chr[1])
        iStr <- as.character(DT[[i]]$strand[1])
        iGA <- geneAnnot[geneAnnot$name == iGene]
        iBlocks <- S4Vectors::mcols(iGA)$blocks %>% unlist %>%
            IRanges::shift(IRanges::start(iGA) - 1)
        tmp <- stringr::str_sub(fastaGenome[[iChr]],
                                start = IRanges::start(iBlocks),
                                end = IRanges::end(iBlocks)) %>%
            paste(collapse = "") %>% Biostrings::DNAString()
        if(iStr == "-"){
            tmp <- Biostrings::reverseComplement(tmp)
        }
        tmp <- stringr::str_split(as.character(tmp), "") %>% unlist
        tibble::add_column(DT[[i]], refSeq = tmp, .after = "txcoor")
    })
}
