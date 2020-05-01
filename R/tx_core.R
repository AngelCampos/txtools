#' txtools: A package to analyze transcriptomic data
#'
#' The txtools package provides functions to analyze genomic data from a
#' transcriptomic perspective. It consists on functions which make the br
#'
#' @section Genome to transcriptome functions:
#' These functions work is to load and process the genomic data into their
#' transcriptomic counterparts.
#' @section Transcriptomic analysis:
#' These functions work is to perform different tasks to analyze the
#' transcriptomic data.
#' @section Meta-analysis:
#' These functions work is to aggretate and analyze transcriptomic data to
#' perform analyses at the meta-transcript level.
#' @section Graphic:
#' Finally graphical functions are available to visualize and inspect final
#' and workflow intermediate results to facilitate the analysis process.
#'
#' @docType package
#' @name txtools
NULL

#' Read paired end bam file by yield size
#'
#' Reads a file in BAM format by blocks of lines equal to a yield size, either
#' automatically calculated or specified by the user, and loads it as a
#' GenomicAlignments object.
#'
#' @param file character. Path to the file to read
#' @param yieldSize numeric. Number of reads to be processed at a time
#' @param scanFlag integer. Flag used to filter reads. See \code{\link[Rsamtools]{ScanBamParam}}
#' @param loadSeq logical. Set to TRUE for loading the sequences contained
#' in the BAM file
#' @param recoverDumpedAligns logical. If set to TRUE ambiguous alignments
#' will be retrieved and output along correctly paired reads. More info:
#' ?GenomicAlignments::readGAlignments
#' @param verbose logical. Set to FALSE to show less information
#'
#' @return GRanges
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
        if(verbose){utils::setTxtProgressBar(pb, i)}
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
#' @author M.A. Garcia-Campos <https://angelcampos.github.io/>
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

#' Load genome
#'
#' Load genome as DNAStrinSet from FASTA file
#'
#' @param fastaFile path to FASTA format file with all genome sequences
#'
#' @return
#' @export
#'
#' @examples
tx_load_genome <- function(fastaFile){
    Biostrings::readDNAStringSet(fastaFile)
}


#' Merge paired-end reads and assign them to gene models
#'
#' @param reads GenomicRanges
#' @param geneAnnot GenomicRanges. Gene annotation loaded via the tx_load_bed()
#' @param overlapType character. Overlap type to filter reads by gene
#' model.
#' @param minReads integer. Minimum number of reads required to overlap a gene
#' model to be part of the output object.
#' @param withSeq logical. Set to TRUE if sequence should be preserved; 'reads'
#' object should contain sequences.
#' @param verbose logical. Set to FALSE to show less information
#'
#' @return
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_reads <- function(reads, geneAnnot, overlapType = "within", minReads = 50,
                     withSeq = F, verbose = T){
    if(class(reads) != "GAlignmentPairs"){
        stop("reads argument should be of class GAlignmentPairs \n")
    }
    if(class(geneAnnot) != "GRanges"){
        stop("geneAnnot argument should be of class GRanges \n")
    }
    if(verbose){
        cat("Processing", length(reads), "reads, using", length(geneAnnot),
            "gene models \n")
    }
    split_i <- hlpr_splitReadsByGenes(reads, geneAnnot, overlapType, minReads)
    geneAnnot <- geneAnnot[which(geneAnnot$name %in% names(split_i))]
    if(length(geneAnnot) > 0){
        if(verbose){
            cat(length(unique(unlist(split_i))), "paired-end reads overlap",
                length(geneAnnot), "gene models \n")
            cat("Filtering reads by gene model... \n")
            if(withSeq){
                cat("Processing sequences. This may take several minutes... \n")
            }
        }
    }else{stop("No genes with overlapped paired-end reads \n")}
    allExons <- exonGRanges(geneAnnot) # All exons in gene models
    OUT <- v_hlpr_ReadsInGene(reads = reads,
                              iGene = geneAnnot$name,
                              geneAnnot = geneAnnot,
                              split_i = split_i,
                              allExons = allExons,
                              withSeq = withSeq,
                              minReads = minReads)
    names(OUT) <- geneAnnot$name
    OUT <- OUT[lapply(OUT, length) %>% unlist %>%
                   magrittr::is_greater_than(minReads)] %>%
        GenomicRanges::GRangesList(compress = TRUE)
    if(verbose){
        cat("Output contains:", lapply(OUT, names) %>% unlist %>% unique %>% length,
            "unique reads in", length(OUT), "gene models \n")
    }
    return(OUT)
}

#' Merge paired-end reads and assign them to gene models (Multicore)
#'
#' @param reads GAlignmentPairs. Paired end genomic alignments to be processed
#' @param geneAnnot GenomicRanges. Gene annotation loaded via the tx_load_bed()
#' @param nCores integer. Number of cores to be used.
#' @param overlapType character. Overlap type to filter reads by gene model.
#' @param minReads integer. Minimum number of reads required to overlap a gene
#' model to be part of the output object.
#' @param withSeq logical. Set to TRUE if sequence should be preserved; 'reads'
#' object should contain sequences.
#' @param verbose logical. Set to FALSE to show less information.
#'
#' @return
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_reads_mc <- function(reads, geneAnnot, nCores, overlapType = "within",
                        minReads = 50, withSeq = F, verbose = T){
    if(class(reads) != "GAlignmentPairs"){
        stop("reads argument should be of class GAlignmentPairs \n")
    }
    if(class(geneAnnot) != "GRanges"){
        stop("geneAnnot argument should be of class GRanges \n")
    }
    if(.Platform$OS.type != "unix"){
        stop("This functions is only available for UNIX operative systems \n")
    }
    if(verbose){
        cat("Processing", length(reads), "paired-end reads, using",
            length(geneAnnot), "gene models \n")
    }
    split_i <- hlpr_splitReadsByGenes(reads, geneAnnot, overlapType, minReads)
    geneAnnot <- geneAnnot[which(geneAnnot$name %in% names(split_i))]
    if(length(geneAnnot) > 0){
        if(verbose){
            cat(length(unique(unlist(split_i))), "paired-end reads overlap",
                length(geneAnnot), "gene models \n")
            cat("Filtering reads by gene model... \n")
            if(withSeq){
                cat("Processing sequences. This may take several minutes... \n")
            }
        }
    }else{
        stop("No genes with overlapped paired-end reads \n")
    }
    allExons <- exonGRanges(geneAnnot)
    OUT <- parallel::mclapply(mc.cores = nCores, geneAnnot$name, function(iGene){
        txtools:::hlpr_ReadsInGene(reads = reads,
                                   iGene = iGene,
                                   geneAnnot = geneAnnot,
                                   split_i = split_i,
                                   allExons = allExons,
                                   withSeq = withSeq,
                                   minReads = minReads)
    })
    names(OUT) <- geneAnnot$name
    OUT <- OUT[lapply(OUT, length) %>% unlist %>% magrittr::is_greater_than(minReads)] %>%
        GenomicRanges::GRangesList(compress = TRUE)
    if(verbose){
        cat("Output contains:", lapply(OUT, names) %>% unlist %>% unique %>% length,
            "unique reads in", length(OUT), "gene models \n")
    }
    return(OUT)
}

#' Filter ranges by a maximum width
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via the tx_reads() or tx_reads_mc()
#' functions.
#' @param thr numeric. Threshold for maximum width size allowed on output.
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
    }) %>% GenomicRanges::GRangesList(compress = TRUE) %>% magrittr::set_names(names(x))
}

#' Calculate coverage table: coverage, 5prime-starts, and 3prime-ends
#'
#' @param x CompressedGRangesList
#'
#' @return data.table
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
            data.table::data.table() %>%
            magrittr::set_rownames(1:len[i]) %>%
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
#' alignments data by gene. Constructed via the tx_reads() or tx_reads_mc()
#' functions.
#' @param nCores integer. Number of cores to use.
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
    cov <- tx_coverage(x) %>% parallel::mclapply(as.vector, mc.cores = nCores)
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
#' alignments data by gene. Constructed via the tx_reads() or tx_reads_mc()
#' functions.
#' @param simplify_IUPAC character. Method to simplify ambiguous reads:
#' 1) 'not': Ambiguous reads will be left as their IUPAC_ambiguous code, e.g.
#' for positions in which a 'G' and and 'A' where read an 'R' will note this.
#' 2) "splitHalf': Ambiguous reads will be divided in half, in odd cases having
#' fractions of reads assigned.
#' 3) "splitForceInt": Ambiguous reads will be divided in half, but forced to be
#' integers, unassigned fraction of reads will be summed and assigned as 'N'.
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
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via the tx_reads() or tx_reads_mc()
#' functions.
#' @param geneAnnot GenomicRanges. Gene annotation loaded via the tx_load_bed()
#'
#' @return
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_genCoorTab <- function(x, geneAnnot){
    if(all(names(x) %in% geneAnnot$name)){
        lapply(names(x), function(iGene){
            tmp2 <- geneAnnot[which(geneAnnot$name == iGene)]
            tmp3 <- c(GenomicAlignments::seqnames(tmp2),
                      GenomicRanges::strand(tmp2)) %>% as.character() %>% c(iGene)
            tmpDT <- rep(tmp3, GenomeInfoDb::seqlengths(x[[iGene]])[iGene]) %>%
                matrix(ncol = 3, byrow = T) %>%
                cbind(exonBlockGen(iGene, geneAnnot)) %>%
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
        stop("Names of x are not contained in geneAnnot$name")
    }
}

#' Summarized Coverage data.table
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via the tx_reads() or tx_reads_mc()
#' functions.
#' @param geneAnnot GenomicRanges. Gene annotation loaded via the tx_load_bed()
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
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via the tx_reads() or tx_reads_mc()
#' functions.
#' @param geneAnnot GenomicRanges. Gene annotation loaded via the tx_load_bed()
#' @param simplify_IUPAC string. Available options are :
#' \itemize{
#' \item "not": Will output the complete nucleotide frequency table including
#' ambiguous reads using the IUPAC ambiguity code.
#' See: \code{\link[Biostrings]{IUPAC_CODE_MAP}}
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
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via the tx_reads() or tx_reads_mc()
#' functions.
#' @param geneAnnot GenomicRanges. Gene annotation loaded via the tx_load_bed()
#' function.
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
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_covNucFreqDT <- function(x, geneAnnot, simplify_IUPAC = "splitForceInt"){
    hlp_cbind3Tabs(tx_genCoorTab(x, geneAnnot),
                   tx_coverageTab(x),
                   tx_nucFreqTab(x, simplify_IUPAC))
}

#' Quantifies reads by gene model from a tx_reads list
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via the tx_reads() or tx_reads_mc()
#' functions.
#'
#' @return table
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_counts <- function(x){
    suppressWarnings(unlist(x)) %>% GenomeInfoDb::seqnames() %>% table()
}

#' Add startRatio to data.table
#'
#' Add read-starts ratio over coverage.
#'
#' @param x data.table. Output of txtools data.tables with coverage information
#' as output of tx_coverageDT() and tx_covNucFreqDT()
#'
#' @return data.table
#' @export
#'
#' @examples
tx_add_startRatio <- function(x){
    tibble::add_column(x, startRatio = x$start_5p / x$cov, .after = "start_5p")
}

#' Add endRatio to data.table
#'
#' Add read-ends ratio over coverage.
#'
#' @param x data.table. Output of txtools data.tables with coverage information
#' as output of tx_coverageDT() and tx_covNucFreqDT()
#'
#' @return data.table
#' @export
#'
#' @examples
tx_add_endRatio <- function(x){
    tibble::add_column(x, endRatio = x$end_3p / x$cov, .after = "end_3p")
}


#' Add site annotation
#'
#' Add a logical variable column in which desired genomic coordinates are
#' marked as TRUE. Useful for example when distinguishing between known
#' modified sites.
#'
#' @param x data.table. Output of txtools data.tables with coverage information
#' as output of tx_coverageDT() and tx_covNucFreqDT()
#' @param GR GenomicRanges. Length 1 ranges which want to be marked in the
#' data.tables objects.
#' @param type character. Type of variable to be added:
#' 1) 'logical': Found coordinates will be marked as TRUE while the rest will be
#' left as FALSE.
#' @param colName character. Name of the new column to be added.
#'
#' @return data.table
#' @export
#'
#' @examples
tx_add_siteAnnotation <- function(x, GR, type = "logical", colName){
    if(class(GR) != "GRanges"){stop("GR must be of class GRanges")}
    if(class(x)[1] != "data.table"){stop("x must be of class data.table")}
    if(!all(GenomicRanges::start(subGR) == GenomicRanges::end(subGR))){
        stop("start and ends are not the same in GR, only 1-nuc-long sites allowed")
    }
    subGR <- GR[as.character(x$chr[1]) == GenomicRanges::seqnames(GR)]
    oNames <- names(x)
    # Logical variable case
    if(type == "logical"){
        addAnnot <- rep(FALSE, nrow(x))
        if(length(subGR) == 0){
            tibble::add_column(x, addAnnot) %>% magrittr::set_names(c(oNames, colName))
        }
        foundGenLoc <- GenomicRanges::start(subGR)[which(GenomicRanges::start(subGR) %in% x$gencoor)]
        if(length(foundGenLoc) == 0){
            tibble::add_column(x, addAnnot) %>% magrittr::set_names(c(oNames, colName))
        }
        addAnnot[match(foundGenLoc, x$gencoor)] <- TRUE
        tibble::add_column(x, addAnnot) %>% magrittr::set_names(c(oNames, colName))
    }
}

#' Add reference sequence to data.table
#'
#' @param DT data.table. A summarized data.table object. See tx_coverageDT(),
#' tx_nucFreq() and tx_covNucFreqDT() functions.
#' @param fastaGenome list. The full reference genome sequences, as prepackaged
#' by BSgenome. See ?BSgenome::available.genomes()
#' @param geneAnnot GRanges. Gene annotation loaded as a GenomicRanges object,
#'  see tx_load_bed().
#'
#' @return data.table
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_add_refSeqDT <- function (DT, fastaGenome, geneAnnot){
    if (class(DT)[1] != "data.table") {
        stop("DT must be of class list")
    }
    iGene <- as.character(unique(DT$gene))
    if (length(iGene) != 1) {
        stop("DT must be a data.table pertaining only one gene")
    }
    iChr <- as.character(DT$chr[1])
    iStr <- as.character(DT$strand[1])
    iGA <- geneAnnot[geneAnnot$name == iGene]
    iBlocks <- S4Vectors::mcols(iGA)$blocks %>% unlist %>%
        IRanges::shift(IRanges::start(iGA) - 1)
    tmp <- stringr::str_sub(fastaGenome[[iChr]], start = IRanges::start(iBlocks),
                            end = IRanges::end(iBlocks)) %>% paste(collapse = "") %>%
        Biostrings::DNAString()
    if (iStr == "-") {
        tmp <- Biostrings::reverseComplement(tmp)
    }
    tmp <- stringr::str_split(as.character(tmp), "") %>%
        unlist
    tibble::add_column(DT, refSeq = tmp, .after = "txcoor")
}

#' Merge lists of data.tables
#'
#' @param DTL1 data.table
#' @param DTL2 data.table
#' @param colsToAdd character. Numeric column(s) to be aggregated (added).
#'
#' @return list
#' @export
#'
#' @examples
tx_aggregate_DTlist <- function(DTL1, DTL2, colsToAdd){
    allNames <- union(names(DTL1), names(DTL2))
    namesInBoth <- intersect(names(DTL1), names(DTL2))
    namesOnly1 <- setdiff(names(DTL1), names(DTL2))
    namesOnly2 <- setdiff(names(DTL2), names(DTL1))
    tmpDTL <- lapply(namesInBoth, function(iGene){
        if(!identical(DTL1[[iGene]][,-..colsToAdd], DTL2[[iGene]][,-..colsToAdd])){
            stop(paste("Coordinate data in", iGene, " data.table is not compatible"))
        }
        tmp <- DTL1[[iGene]][,..colsToAdd] + DTL2[[iGene]][,..colsToAdd]
        cbind(DTL1[[iGene]][,-..colsToAdd], tmp)
    }) %>% magrittr::set_names(namesInBoth)
    c(tmpDTL, DTL1[namesOnly1], DTL2[namesOnly2])[allNames]
}


#' Merge data.tables in list to a single data.table
#'
#' @param x list. List of data.tables. A summarized data.table object.
#' See tx_coverageDT(), tx_nucFreqDT() and tx_covNucFreqDT() functions.
#'
#' @return data.table
#' @export
#'
#' @examples
tx_merge_DT <- function(x){
    do.call(x, what = rbind)
}


#' Split data.table to list of data.tables
#'
#' Split data.table back to list with individual data.tables by 'gene' names
#'
#' @param x data.table. Merged data.table as output by tx_merge_DT()
#'
#' @return list
#' @export
#'
#' @examples
tx_split_DT <- function(x){
    tmp <- split(x, by = "gene")
    lapply(tmp, function(y){
        y[order(y$txcoor),]
    })
}
