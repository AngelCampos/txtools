#' txtools: A package to analyze transcriptomic data
#'
#' The txtools package provides functions to analyze genomic data from a
#' transcriptomic perspective. It consists on functions which make the br
#'
#' @section tx_load_*():
#' These functions work is to load and process the genomic data into their
#' transcriptomic counterparts.
#' @section tx_reads():
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
#' @aliases txtools-package
#' @keywords DataImport DataRepresentation Coverage Epitranscriptomics RNASeq
#' @keywords Transcription SNP
"_PACKAGE"
NULL

# Loading files into R #########################################################

#' Read paired end bam file by yield size
#'
#' Reads a file in BAM format by blocks of lines equal to a yield size, either
#' automatically calculated or specified by the user, and loads it as a
#' GenomicAlignments object.
#'
#' @param file character. Path to the file to read
#' @param pairedEnd logical. Set to FALSE if reads in BAM file are single-end,
#' set to TRUE if reads are paired-end.
#' @param yieldSize numeric. Number of reads to be processed at a time
#' @param scanFlag integer. Flag used to filter reads.
#' See \code{\link[Rsamtools]{ScanBamParam}}
#' @param loadSeq logical. Set to TRUE for loading the sequences contained
#' in the BAM file
#' @param verbose logical. Set to FALSE to show less information
#'
#' @return GRanges
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
#' # Loading in-package BAM file
#' bamFile <- system.file("extdata", "example_hg19.bam", package = "txtools")
#' hg19_bam <- tx_load_bam(bamFile, pairedEnd = TRUE, loadSeq = TRUE, verbose = TRUE)
#' summary(hg19_bam)
tx_load_bam <- function(file, pairedEnd, yieldSize = 100000,
                        scanFlag = "default", loadSeq = FALSE, verbose = TRUE){
    if(!is.logical(pairedEnd)){stop("Argument 'pairedEnd' must be of class logical")}
    if(!is.logical(loadSeq)){stop("Argument 'loadSeq' must be of class logical")}
    if(verbose){cat("Reading number of records in file \n")}
    #Remove optical duplicates, low quality reads, low quality reads, non-paired
    if(scanFlag == "default" & pairedEnd){
        scanFlag <- Rsamtools::scanBamFlag(isDuplicate = FALSE,
                                           isNotPassingQualityControls = FALSE,
                                           isPaired = TRUE)
    }else if(scanFlag == "default" & !pairedEnd){
        scanFlag <- Rsamtools::scanBamFlag(isDuplicate = FALSE,
                                           isNotPassingQualityControls = FALSE)
    }
    if(loadSeq){
        readParams <- Rsamtools::ScanBamParam(flag = scanFlag, what = "seq")
    }else{
        readParams <- Rsamtools::ScanBamParam(flag = scanFlag)
    }
    # CountBam records (each paired end alignments equals two records)
    bC <- Rsamtools::countBam(file)
    if(bC$records == 0){stop("BAM file is empty \n")}
    if(verbose){cat(bC$records, "number of BAM records \n")}
    if(verbose){
        cat("Loading BAM file \n")
        pb <- utils::txtProgressBar(style = 3)
    }
    BAMFILE <- Rsamtools::BamFile(file, yieldSize = yieldSize)
    readCycles <- 1:ceiling(bC$records/(yieldSize * ifelse(pairedEnd, 2, 1)))
    pbJumps <- seq(0, 1, by = min(1/(length(readCycles)-1), 1))
    #  Opens and reads BAM files in chunks of length = yieldSize
    open(BAMFILE)
    reads <- lapply(pbJumps, function(i){
        if(verbose){utils::setTxtProgressBar(pb, i)}
        if(pairedEnd){
            suppressWarnings(
                GenomicAlignments::readGAlignmentPairs(BAMFILE,
                                                       use.names = TRUE,
                                                       param = readParams))
        }else if(pairedEnd == FALSE){
            suppressWarnings(
                GenomicAlignments::readGAlignments(BAMFILE,
                                                   use.names = TRUE,
                                                   param = readParams))
        }
    }) %>% do.call(what = c)
    close(BAMFILE)
    if(verbose){close(pb)}
    if(loadSeq){reads <- hlp_cleanBam_emptySeq(reads, verbose)}
    if(verbose){cat(" \n")}
    if(verbose){
        cat(length(reads), "Reads succesfully loaded \n")
        cat("Dumped ambiguous reads:", GenomicAlignments::countDumpedAlignments(), "\n")
    }
    if(GenomicAlignments::countDumpedAlignments() > 0){
        warning(GenomicAlignments::countDumpedAlignments(), " dumped ambiguous reads. ",
                "Use 'getDumpedAlignments()' to retrieve ",
                "them from the dump environment.")
    }
    return(reads)
}

#' Load gene models from bed-12 and bed-6 files
#'
#' Reads and loads bed files checking for compliance with bed file structure
#'
#' @param bedfile character. Gene annotation file in bed-6 or bed-12 format
#'
#' @return GRanges
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
#' Load genome as DNAStrinSet from FASTA file. Alias for
#' \code{\link[Biostrings]{readDNAStringSet}}
#'
#' @param fastaFile path to FASTA format file with all genome sequences
#'
#' @return DNAStringSet
#' @export
#'
#' @examples
tx_load_genome <- function(fastaFile){
    Biostrings::readDNAStringSet(fastaFile)
}

#' Loading RDS files into data.tables
#'
#' @param file File name or path to .rds file to be loaded
#'
#' @return data.table
#' @export
#'
#' @examples
tx_load_rdsDT <- function(file){
    data.table::data.table(readRDS(file))
}

# Assignment of GenomicAlignments to genes #####################################

#' Transcriptomic reads convertion
#'
#' Assign aligned reads to their respective gene models and convert their
#' positions into their corresponding transcriptomic coordinate system.
#' It also stitches together paired-end aligned reads into a single sequence
#' in which dots '.' separate Read1 and Read2 by their corresponding insert.
#'
#' To retrieve unassigned alignments use the function tx_getUnassignedAlignments()
#'
#' @param reads GAlignments or GAlignmentPairs. Genomic alignments to be processed
#' @param geneAnnot GenomicRanges. Gene annotation loaded via the tx_load_bed()
#' @param minReads integer. Minimum number of reads required to overlap a gene
#' @param withSeq logical. Set to TRUE if sequence should be preserved; 'reads'
#' object should contain sequences.
#' @param verbose logical. Set to FALSE to show less information.
#' @param nCores integer. Number of cores to use to run function. Multi-core
#' capability not available in Windows OS.
#'
#' @aliases tx_reads_mc
#' @aliases tx_flushUnassigned
#'
#' @return GRanges
#' @export
tx_reads <- function(reads, geneAnnot, minReads = 50, withSeq = FALSE, verbose = TRUE,
                     nCores = 1){
    # Checks
    tx_flushUnassigned()
    check_mc_windows(nCores)
    check_integerGreaterThanZero_arg(minReads, "minReads")
    check_GA_reads_compatibility(reads, geneAnnot)
    if(withSeq){check_BAM_has_seq(reads)}
    if(!class(reads) %in% c("GAlignmentPairs", "GAlignments")){
        stop("reads argument should be of class GAlignmentPairs. \n")
    }
    if(class(geneAnnot) != "GRanges"){
        stop("geneAnnot argument should be of class GRanges. \n")
    }
    pairedEnd <- switch(class(reads), GAlignmentPairs = TRUE, GAlignments = FALSE)
    # Main program
    if(verbose){
        cat("Processing", length(reads), "reads, using", length(geneAnnot),
            "gene models. \n")
    }
    # Spliting reads by gene
    overlapType <-  "within"
    split_i <- switch(class(reads),
                      GAlignmentPairs = hlpr_splitReadsByGenes(reads, geneAnnot,
                                                               overlapType, minReads),
                      GAlignments = hlpr_splitReadsByGenes_singleEnd(reads, geneAnnot,
                                                                     overlapType, minReads))
    geneAnnot <- geneAnnot[which(geneAnnot$name %in% names(split_i))]
    if(length(geneAnnot) > 0){
        if(verbose){
            cat(length(unique(unlist(split_i))), "alignments overlap",
                length(geneAnnot), "gene models \n")
            cat("Assigning alignments to gene model... \n")
            if(withSeq){
                cat("Processing sequences. This may take several minutes",
                    "depending on geneAnnot size ... \n")
            }
        }
    }else{stop("There where no reads, or not enough reads overlapping the gene models. \n")}
    allExons <- exonGRanges(geneAnnot) # All exons in gene models
    # Pass from genomic coordinates into transcriptomic coordinates and
    # Stitch sequences for paired-end reads
    if(pairedEnd){
        OUT <- parallel::mclapply(mc.cores = nCores, geneAnnot$name, function(iGene){
            hlpr_ReadsInGene(reads = reads,
                             iGene = iGene,
                             geneAnnot = geneAnnot,
                             split_i = split_i,
                             allExons = allExons,
                             withSeq = withSeq,
                             minReads = minReads)
        })
    }else if(!pairedEnd){
        OUT <- parallel::mclapply(mc.cores = nCores, geneAnnot$name, function(iGene){
            hlpr_ReadsInGene_SingleEnd(reads = reads,
                                       iGene = iGene,
                                       geneAnnot = geneAnnot,
                                       split_i = split_i,
                                       allExons = allExons,
                                       withSeq = withSeq,
                                       minReads = minReads)
        })
    }
    names(OUT) <- geneAnnot$name
    OUT <- OUT[lapply(OUT, length) %>% unlist %>% magrittr::is_greater_than(minReads)] %>%
        GenomicRanges::GRangesList()
    if(verbose){
        cat("Output contains:", lapply(OUT, names) %>% unlist %>% unique %>% length,
            "unique reads in", length(OUT), "gene models \n")
    }
    # Dump not assigned alignments into dump environment
    if(length(setdiff(names(reads), unique(unlist(lapply(OUT, names))))) > 1L){
        warning("Some alignments were not assigned to any gene, you can ",
                "retrieve them using the tx_getUnassignedAlignments() function.")
        hlp_dump_notAssigned(reads, OUT)
    }
    return(OUT)
}

#' Retrieve dumped alignments
#'
#' Retrieve the read mappings that were not assigned to any gene by
#' \code{\link{tx_reads}}
#'
#' @return
#' @export
#'
tx_getUnassignedAlignments <- function(){
    objnames <- ls(envir = .dumpEnvirTxtools())
    do.call(c, unname(mget(objnames, envir = .dumpEnvirTxtools())))
}

# Manipulating GenomicRanges ###################################################

#' Extending GRanges 5' and 3' UTR blocks
#'
#' A function that extends the 5' and/or 3' UTR regions of GRanges including the
#' 'blocks' column, which defines the exon structure of a transcript.
#' As those loaded by ´tx_load_bed()´.
#'
#' @param GR GRanges GRanges containing a blocks
#' @param ext_5p integer Number of bp for the 5' UTR blocks to be extended
#' @param ext_3p integer Number of bp for the 3' UTR blocks to be extended
#'
#' @return GRanges
#' @export
#'
#' @examples
tx_extend_UTR <- function(GR, ext_5p = 0, ext_3p = 0){
    if(is.null(GR$blocks)){
        stop("GR does not containg 'blocks' meta column.")
    }
    GR_out <- plyranges::stretch(plyranges::anchor_5p(GR), ext_3p)
    GR_out <- plyranges::stretch(plyranges::anchor_3p(GR_out), ext_5p)
    GR_out$blocks <- stretchBlocks_3p(GR_out$blocks, ext_3p, GenomicRanges::strand(GR))
    GR_out$blocks <- stretchBlocks_5p(GR_out$blocks, ext_5p, GenomicRanges::strand(GR))
    return(GR_out)
}

# Manipulate GRangesList #######################################################

#' Filter ranges by a maximum width
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via tx_reads().
#' @param thr integer. Threshold for maximum width size allowed on output.
#' @param nCores integer. Number of cores to use to run function. Multicore
#' capability is not available in Windows OS.
#'
#' @return CompressedGRangesList
#' @export
#'
#' @author M.A. Garcia-Campos
#' @aliases tx_filter_max_width
#' @examples
tx_filter_maxWidth <- function(x, thr, nCores = 1){
    check_integerGreaterThanZero_arg(nCores, "nCores")
    check_integer_arg(thr, "thr")
    check_mc_windows(nCores)
    tmp <- GenomicAlignments::width(x) %>% magrittr::is_weakly_less_than(thr)
    parallel::mclapply(mc.cores = nCores, seq(1, length(x)), function(i){
        if(all(tmp[[i]])){
            x[[i]]
        }else{
            x[[i]][tmp[[i]]]
        }
    }) %>% GenomicRanges::GRangesList() %>% magrittr::set_names(names(x))
}


#' Sampling alignments
#'
#' Sampling of alignmnets in a GRanges list using a binomial distribution.
#'
#' @param x CompressedGRangesList. A list containing alignments, meant to be used
#' for the output of the \code{\link{tx_reads}} function.
#' @param p Probabilty for each read to be sampled.
#' @param nCores integer. Number of cores to use to run function. Multicore
#' capability not available in Windows OS.
#'
#' @return CompressedGRangesList.
#' @export
#'
#' @examples
tx_sample_GRList <- function(x, p, nCores = 1){
    parallel::mclapply(mc.cores = nCores, x, function(GR){
        GR[sample(c(T,F), prob = c(p, 1-p), replace = T, size = length(GR))]
    }) %>% GenomicRanges::GRangesList()
}


# Generating txtools DTs #######################################################

#' Summarized Coverage data.table
#'
#' This function constructs a data.table that contains coverage metrics
#' per nucleotide for all transcripts with corresponding GRanges in its first argument 'x':
#'\itemize{
#'\item cov = Insert coverage
#'\item start_5p = read-start counts
#'\item end_3p = read-end counts
#'}
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via the \code{\link{tx_reads}}
#' function.
#' @param geneAnnot GenomicRanges. Gene annotation loaded via the \code{\link{tx_load_bed}}
#' @param genome list. The full reference genome sequences, as prepackaged
#' by BSgenome, See ?BSgenome::available.genomes(); or loaded by \code{\link{tx_load_genome}}
#' @param fullDT logical. Set to TRUE if it is desired to output a data.table
#' with all genes and in the same order as 'geneAnnot' object.
#' @param nCores integer. Number of cores to use to run function.
#'
#' @return data.table
#' @export
#'
#' @family makeDT functions
#' @author M.A. Garcia-Campos
#' @aliases tx_coverageDT
#'
#' @examples
tx_makeDT_coverage <- function(x, geneAnnot, genome = NULL, fullDT = FALSE,
                               nCores = 1){
    tx_flushUnassigned()
    check_mc_windows(nCores)
    check_integerGreaterThanZero_arg(nCores, "nCores")
    OUT <- hlp_cbind2Tabs(hlp_genCoorTab_mc(x, geneAnnot, nCores),
                          hlp_coverageTab_mc(x, nCores)) %>% tx_merge_DT()
    if(fullDT){
        OUT <- tx_complete_DT(DT = OUT, geneAnnot = geneAnnot, genome = genome, nCores = nCores)
    }
    if(!is.null(genome)){
        tx_add_refSeqDT(OUT, genome = genome, geneAnnot = geneAnnot, nCores = nCores)
    }else{return(OUT)}
}

#' Summarized Nucleotide Frequency data.table
#'
#' This function constructs a list of data.tables that contains nucleotide frequency
#' metrics per nucleotide by transcript:
#'\itemize{
#'\item A = Adenine
#'\item C = Cytosine
#'\item G = Guanine
#'\item T = Thymine
#'\item N = Undetermined nucleotide
#'\item - = Deletion
#'\item . = Insert, not read gap between read1 and read2
#'}
#' The function requires the input of a GRangesList object output by the
#' \code{\link{tx_reads}} function, which should contain sequence alignments in the
#' transcriptomic space, and a gene annotation in GRanges format, as loaded by
#' the \code{\link{tx_load_bed}} function.
#'
#' This function allows for usage of multiple cores to reduce processing times
#' in UNIX-like OS.
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via the \code{\link{tx_reads}} function.
#' @param geneAnnot GenomicRanges. Gene annotation loaded via the tx_load_bed()
#' @param genome list. The full reference genome sequences, as prepackaged
#' by BSgenome, See ?BSgenome::available.genomes(); or loaded by \code{\link{tx_load_genome}}
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
#' @param fullDT logical. Set to TRUE if it is desired to output a data.table
#' with all genes and in the same order as 'geneAnnot' object.
#' @param nCores integer. Number of cores to use to run function.
#'
#' @return data.table
#' @export
#'
#' @family makeDT functions
#' @author M.A. Garcia-Campos
#' @aliases tx_nucFreqDT
#'
#' @examples
tx_makeDT_nucFreq <- function(x, geneAnnot, genome = NULL,
                              simplify_IUPAC = "splitForceInt", fullDT = FALSE,
                              nCores = 1){
    tx_flushUnassigned()
    check_mc_windows(nCores)
    check_integerGreaterThanZero_arg(nCores, "nCores")
    check_GR_has_seq(x, "x")
    OUT <- hlp_cbind2Tabs(hlp_genCoorTab_mc(x, geneAnnot, nCores),
                          hlp_nucFreqTab_mc(x, simplify_IUPAC, nCores)) %>% tx_merge_DT()
    if(fullDT){
        OUT <- tx_complete_DT(DT = OUT, geneAnnot = geneAnnot, genome = genome, nCores = nCores)
    }
    if(!is.null(genome)){
        tx_add_refSeqDT(OUT, genome = genome, geneAnnot = geneAnnot, nCores = nCores)
    }else{return(OUT)}
}

#' Summarized Coverage & Nucleotide Frequency data.table
#'
#' This function constructs a list of data.tables that contains nucleotide frequency
#' metrics per nucleotide by transcript:
#'\itemize{
#'\item REF = Equal to reference genome
#'\item . = Insert. Not read gap between read1 and read2, counted for coverage
#' computation
#'\item A = Adenine
#'\item C = Cytosine
#'\item G = Guanine
#'\item T = Thymine
#'\item - = Deletion
#'\item N = Undetermined nucleotide
#'}
#' The function requires the input of a GRangesList object output by the
#' \code{\link{tx_reads}} function, which should contain sequence alignments in the
#' transcriptomic space, and a gene annotation in GRanges format, as loaded by
#' the \code{\link{tx_load_bed}} function.
#'
#' This function allows for usage of multiple cores to reduce processing times
#' in UNIX-like OS.
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via tx_reads().
#' @param geneAnnot GenomicRanges. Gene annotation loaded via the tx_load_bed()
#' function.
#' @param genome list. The full reference genome sequences, as prepackaged
#' by BSgenome, See ?BSgenome::available.genomes(); or loaded by \code{\link{tx_load_genome}}
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
#' @param fullDT logical. Set to TRUE if it is desired to output a data.table
#' with all genes and in the same order as 'geneAnnot' object.
#' @param nCores integer. Number of cores to use to run function.
#'
#' @return data.table
#' @export
#'
#' @family makeDT functions
#' @author M.A. Garcia-Campos
#' @aliases tx_covNucFreqDT
#'
#' @examples
tx_makeDT_covNucFreq <- function(x, geneAnnot, genome = NULL,
                                 simplify_IUPAC = "splitForceInt", fullDT = FALSE,
                                 nCores = 1){
    tx_flushUnassigned()
    check_GR_has_seq(x, "x")
    check_mc_windows(nCores)
    check_integerGreaterThanZero_arg(nCores, "nCores")
    OUT <- hlp_cbind3Tabs(hlp_genCoorTab_mc(x, geneAnnot, nCores),
                          hlp_coverageTab_mc(x, nCores),
                          hlp_nucFreqTab_mc(x, simplify_IUPAC, nCores)) %>% tx_merge_DT()
    if(fullDT){
        OUT <- tx_complete_DT(DT = OUT, geneAnnot = geneAnnot, genome = genome, nCores = nCores)
    }
    if(!is.null(genome)){
        tx_add_refSeqDT(OUT, genome = genome, geneAnnot = geneAnnot, nCores = nCores)
    }else{return(OUT)}
}

# Manipulating data.tables and DT lists ########################################

#' Merge data.tables in list to a single data.table
#'
#' @param DTL list. A data.tables list. See tx_coverageDT(), tx_nucFreqDT()
#' and tx_covNucFreqDT() functions.
#'
#' @return data.table
#' @export
#'
tx_merge_DT <- function(DTL){
    do.call(DTL, what = rbind)
}


#' Split data.table to list of data.tables
#'
#' Split data.table back to list with individual data.tables by 'gene' names
#'
#' @param DT data.table. Merged data.table as output by tx_merge_DT()
#' @param dropEmpty logical. Drops empty list elements, which occurs when data
#' of genes have been entirely removed, but kept listed in the x$gene factor levels.
#'
#' @return list
#' @export
#'
tx_split_DT <- function(DT, dropEmpty = TRUE){
    tmp <- split(DT, f = DT$gene, drop = dropEmpty)
    # lapply(tmp, function(y) {
    #     y[order(y$txcoor), ]
    # })
    tmp
}

#' Cutting 5' and 3' ends of data.table using txcoors
#'
#' This function removes the heading (5'UTR) and trailing (3' UTR) portions of
#' data.tables.
#'
#' @param DT data.table A data.table as generated by the \code{\link{tx_coverageDT}},
#' \code{\link{tx_nucFreqDT}}, and \code{\link{tx_covNucFreqDT}} functions.
#' @param cut_5p integer Basepairs to me removed from the start of DT for each gene
#' @param cut_3p integer Basepairs to me removed from the end of DT for each gene
#'
#' @return data.table
#' @export
#'
#' @examples
tx_cutEnds_DT <- function(DT, cut_5p = 0, cut_3p = 0){
    DTL <- tx_split_DT(check_DT(DT))
    lapply(DTL, function(x){
        hlp_remove_UTR(x, cut_5p, cut_3p)
    }) %>% tx_merge_DT()
}

#' Order txDT
#'
#' @param DT
#'
#' @return
#' @export
#'
#' @examples
tx_orderDT <- function(DT){
    DT[order(DT$gene, DT$txcoor),]
}

#' Complete a DT object missing genes
#'
#' A function that adds the genomic and transcriptomic coordinates of genes
#' with no data in a data.table generated by the \code{\link{tx_coverageDT}},
#' \code{\link{tx_nucFreqDT}}, and \code{\link{tx_covNucFreqDT}} functions.
#'
#' This function is meant to equalize datasets which don't contain some genes,
#' as such the order of the genes is the same as in the GRanges object used as
#' gene annotation.
#'
#' @param DT data.table
#' @param geneAnnot GRanges. Gene annotation, gene models, as loaded using the
#' \code{\link{tx_load_bed}} function.
#' @param genome list. The full reference genome sequences, as loaded by
#' \code{\link{tx_load_genome}}; or prepackaged by BSgenome, see
#' \code{\link[BSgenome]{available.genomes}}
#' @param nCores numeric
#'
#' @return data.table
#' @export
#'
#' @examples
tx_complete_DT <- function(DT, geneAnnot, genome = NULL, nCores = 1){
    check_mc_windows(nCores)
    DT <- check_DT(DT)
    if(!all(unique(DT$gene) %in% geneAnnot$name)){
        stop("All genes in DT must be in geneAnnot, missing genes from ",
             "geneAnnot will be added")
    }
    missGenes <- setdiff(geneAnnot$name, unique(DT$gene))
    if(length(missGenes) == 0){
        warning("No missing genes, tx_complete_DT() didn't add genes.")
        return(DT)
    } # No missing genes
    tmpCoorTabs <- hlpr_genCoorTabGenes(genes = missGenes,
                                        geneAnnot = geneAnnot,
                                        nCores = nCores) %>%
        tx_merge_DT()
    # Add refSeq if present in DT
    if("refSeq" %in% names(DT)){
        tmpCoorTabs <- tx_add_refSeqDT(DT = tmpCoorTabs, genome = genome,
                                       geneAnnot = geneAnnot, nCores = nCores)
    }
    missCols <- setdiff(names(DT), names(tmpCoorTabs))
    tmpCoorTabs <- cbind(tmpCoorTabs,
                         data.table::data.table(matrix(0,
                                                       nrow = nrow(tmpCoorTabs),
                                                       ncol = length(missCols))) %>%
                             magrittr::set_names(missCols)) %>% tx_split_DT()
    DTL <- tx_split_DT(DT)
    completeDTL <- c(DTL, tmpCoorTabs)
    completeDTL <- completeDTL[geneAnnot$name]
    tx_merge_DT(tx_orderDT(completeDTL))
}


#' shift column in txDT
#'
#' @param DT
#' @param colToShift
#' @param direction
#' @param bp
#' @param nCores
#'
#' @return
#' @export
#'
#' @examples
tx_shift_geneWise <- function(DT, colToShift, direction, bp, nCores = 1){
    if(!colToShift %in% colnames(DT)){
        stop("colToShift is not a column of DT")
    }
    DTL <- tx_split_DT(DT)
    if(direction == "upstream"){
        OUT <- parallel::mclapply(mc.cores = nCores, DTL, function(x){
            x[[colToShift]] <- c(utils::tail(x[[colToShift]], -bp), rep(NA, bp))
            return(x)
        }) %>% tx_merge_DT()
    }else if(direction == "downstream"){
        OUT <- parallel::mclapply(mc.cores = nCores, DTL, function(x){
            x[[colToShift]] <- c(rep(NA, bp), utils::head(x[[colToShift]], -bp))
            return(x)
        }) %>% tx_merge_DT()
    }else{
        stop("Input to argument 'direction' is not either 'downstream' or 'upstream'.")
    }
    return(OUT)
}

#' Unify lists of txDTs
#'
#' @param txDTL
#' @param geneAnnot
#' @param genome
#' @param type
#' @param nCores
#'
#' @return
#' @export
#'
#' @examples
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

# data.table functions #########################################################

#' Add reference sequence to a data.table
#'
#' Adds a column of characters representing the corresponsing nucleotide at
#' the reference genome position, considering the strand of the transcript.
#' The input DT must be the output, or resemble it, of the
#'
#'
#' @param DT data.table. A data.table object. See
#' \code{\link{tx_makeDT_coverage}}, \code{\link{tx_makeDT_nucFreq}} or
#' \code{\link{tx_makeDT_covNucFreq}} functions.
#' @param genome list. The full reference genome sequences, as prepackaged
#' by BSgenome, See ?BSgenome::available.genomes(); or loaded by \code{\link{tx_load_genome}}
#' @param geneAnnot GRanges. Gene annotation loaded as a GenomicRanges object,
#'  see tx_load_bed().
#' @param nCores integer. Number of cores to use to run function. Multicore
#' capability not available in Windows OS.
#'
#' @return data.table
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_add_refSeqDT <- function(DT, genome, geneAnnot, nCores = 1){
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("refSeq")
    check_mc_windows(nCores)
    check_integerGreaterThanZero_arg(nCores, "nCores")
    DT <- check_DT(DT)
    check_GA_genome_chrCompat(geneAnnot = geneAnnot, genome = genome)
    DTL <- tx_split_DT(DT)
    tmp <- parallel::mclapply(mc.cores = nCores, DTL, function(DT){
        hlp_add_refSeqDT(DT, genome = genome, geneAnnot = geneAnnot)
    })
    OUT <- tx_merge_DT(tmp)
    return(OUT)
}

#' Add counts of nucleotide reads different to the reference sequence
#'
#' Add a column to a txDT of the sum of nucleotide frequency different to the
#' reference sequence, counting deletions, without considering 'N's nor inserts
#' '.' into the calculation.
#'
#' @param DT data.table. A data.table object. See
#' \code{\link{tx_makeDT_coverage}}, \code{\link{tx_makeDT_nucFreq}} or
#' \code{\link{tx_makeDT_covNucFreq}} functions.
#'
#' @return data.table
#' @export
tx_add_misincCount <- function(DT){
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("misincCount")
    selNucs <- setdiff(intersect(txtools::IUPAC_code_2nucs, names(DT)), c(".", "N"))
    tmp <- DT[, selNucs, with = FALSE]
    OUT <- rep(NA, nrow(DT))
    nucsInRef <- unique(DT$refSeq)
    for(i in nucsInRef){
        selNucs <- setdiff(c("A", "C", "G", "T", "-"), i)
        OUT[which(DT$refSeq == i)] <- rowSums(tmp[which(DT$refSeq == i), selNucs, with = F])
    }
    tibble::add_column(DT, misincCount = OUT)
}


#' Add total number of nucleotide reads
#'
#' Add a column to DT of the sum of nucleotide frequencies for meaningful
#' nucleotide reads not counting undetermined 'N' and inserts '.'.
#'
#' @param DT data.table
#'
#' @return data.table
#' @export
#'
#'
#' @examples
tx_add_nucTotal <- function(DT){
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("nucTotal")
    selNucs <- setdiff(intersect(txtools::IUPAC_code_2nucs, names(DT)), c(".", "N"))
    out <- rowSums(DT[, selNucs, with = FALSE])
    tibble::add_column(DT, nucTotal = out)
}

#' Add misincorporation to total nucleotide reads ratio
#'
#' Add a column to txDT of the ratio of different nucleotides to the total of
#' nucleotide reads, not counting undetermined reads 'N' and inserts '.'.
#'
#' @param DT data.table
#' @param addMisinandTotalCols Set to TRUE to add counts of total nucleotides
#' read (nucTotal) and different to reference nucleotide (misincCount) columns.
#'
#' @return data.table
#' @export
#'
#' @seealso \code{\link{tx_add_diffNucToRef}} and \code{\link{tx_add_nucTotal}}
tx_add_misincRate <- function(DT, addMisinandTotalCols = FALSE){
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("misincRate")
    tmpDT <- tx_add_misincCount(DT) %>% tx_add_nucTotal()
    tmp <- round(tmpDT$misincCount / tmpDT$nucTotal, 6)
    if(addMisinandTotalCols){DT <- tmpDT}
    tibble::add_column(DT, misincRate = tmp)
}
#' Nucleotide specific misincorporation rate
#'
#' Calculates the misincorporation rate from a \bold{ref}erence nucleotide to an
#' specific \bold{mis}incorporated nucleotide. For example, when looking for
#' positions in cytidines that were read as thymine in RNA after bisulphite treatment.
#'
#' @param DT
#' @param minCov
#' @param refNuc Reference nucleotide
#' @param misNuc Misincorporated nucleotide
#'
#' @return
#' @export
#'
#' @examples
tx_add_misincorpRateNucSpec <- function(DT, minCov = 50, refNuc, misNuc){
    DT <- data.table::data.table(DT)
    newColName <- paste0("MR_", refNuc, "to", misNuc)
    if(newColName %in% colnames(DT)){
        DT <- DT[,!grepl(newColName, colnames(DT))]
    }
    tmp <- round(DT[[misNuc]] / (DT[[misNuc]] + DT[[refNuc]]), 6)
    tmp[DT$cov < minCov] <- NA
    tmp[!(DT$refSeq %in% refNuc)] <- NA
    DT[, newColName] <- tmp
    DT
}

#' Add starts to coverage ratio
#'
#' Add a column to DT of the read-starts to coverage ratio.
#'
#' @param DT data.table. Output of txtools data.tables with coverage information
#' as output of \code{\link{tx_coverageDT}} \code{\link{tx_covNucFreqDT}}
#' functions
#' @param minCov numeric. Minimum coverage required to output ratio. If coverage
#' is less then an NA is output in that position.
#'
#' @return data.table
#' @export
#'
tx_add_startRatio <- function(DT, minCov = 50){
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("startRatio")
    tmp <- (DT$start_5p) / (DT$cov)
    tmp[DT$cov < minCov] <- NA
    tibble::add_column(DT, startRatio = tmp)
}

#' Add starts to coverage ratio 1 bp downstream
#'
#' Add a column to DT of the read-starts to coverage ratio, shifted 1 base-pair
#' downstream. This means that the last measurement in any gene is always an NA
#' to account that there was no measurement to be shifted.
#'
#'
#' @param DT data.table. Output of txtools data.tables with coverage information
#' as output of \code{\link{tx_coverageDT}} \code{\link{tx_covNucFreqDT}}
#' functions
#' @param minCov numeric. Minimum coverage required to output ratio. If coverage
#' is less then an NA is output in that position.
#'
#' @return data.table
#' @export
#'
#' @examples
tx_add_startRatio1bpDS <- function(DT, minCov = 50){
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("startRatio1bpDS")
    tmp <- (DT$start_5p) / (DT$cov)
    tmp[DT$cov < minCov] <- NA
    DTL <- tibble::add_column(DT, startRatio1bpDS = tmp) %>% txtools::tx_split_DT()
    DT <- lapply(DTL, function(DT){
        DT$startRatio1bpDS <- c(utils::tail(DT$startRatio1bpDS, -1), NA)
        DT
    }) %>% txtools::tx_merge_DT()
    return(DT)
}

#' Add starts to coverage ratio 1 bp upstream
#'
#' @param DT data.table. Output of txtools data.tables with coverage information
#' as output of \code{\link{tx_coverageDT}} \code{\link{tx_covNucFreqDT}}
#' functions
#' @param minCov numeric. Minimum coverage required to output ratio. If coverage
#' is less then an NA is output in that position.
#'
#' @return data.table
#' @export
#'
#' @examples
tx_add_startRatio1bpUS <- function(DT, minCov = 50){
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("startRatio1bpUS")
    tmp <- (DT$start_5p) / (DT$cov)
    tmp[DT$cov < minCov] <- NA
    DTL <- tibble::add_column(DT, startRatio1bpUS = tmp) %>% txtools::tx_split_DT()
    DT <- lapply(DTL, function(DT){
        DT$startRatio1bpUS <- c(NA, utils::head(DT$startRatio1bpUS, -1))
        DT
    }) %>% txtools::tx_merge_DT()
    return(DT)
}



#' Add ends to coverage ratio 1 bp downstream
#'
#' Adds a column to DT of the read-ends to coverage ratio.
#'
#' @param DT data.table. Output of txtools data.tables with coverage information
#' as output of \code{\link{tx_coverageDT}} \code{\link{tx_covNucFreqDT}}
#' functions
#' @param minCov numeric. Minimum coverage required to output ratio. If coverage
#' is less then an NA is output in that position.
#'
#' @return data.table
#' @export
#'
tx_add_endRatio <- function(DT, minCov = 50){
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("endRatio")
    tmp <- (DT$end_3p) / (DT$cov)
    tmp[DT$cov < minCov] <- NA
    tibble::add_column(DT, endRatio = tmp)
}

#' Add ends to coverage ratio 1bp down-stream
#'
#' Add a column to DT of the read-ends to coverage ratio, shifted 1 base-pair
#' downstream. This means that the last measurement in any gene is always an NA
#' to account that there was no measurement to be shifted.
#'
#' @param DT data.table. Output of txtools data.tables with coverage information
#' as output of \code{\link{tx_coverageDT}} \code{\link{tx_covNucFreqDT}}
#' functions
#' @param minCov numeric. Minimum coverage required to output ratio. If coverage
#' is less then an NA is output in that position.
#'
#' @return data.table
#' @export
#'
tx_add_endRatio1bpDS <- function(DT, minCov = 50){
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("endRatio1bpDS")
    tmp <- (DT$end_3p) / (DT$cov)
    tmp[DT$cov < minCov] <- NA
    DTL <- tibble::add_column(DT, endRatio1bpDS = tmp) %>% txtools::tx_split_DT()
    DT <- lapply(DTL, function(DT){
        DT$endRatio1bpDS <- c(utils::tail(DT$endRatio1bpDS, -1), NA)
        DT
    }) %>% txtools::tx_merge_DT()
}


#' Add ends to coverage ratio 1bp up-stream
#'
#' @param DT data.table. Output of txtools data.tables with coverage information
#' as output of \code{\link{tx_coverageDT}} \code{\link{tx_covNucFreqDT}}
#' functions
#' @param minCov numeric. Minimum coverage required to output ratio. If coverage
#' is less then an NA is output in that position. To output all positions assign 0.
#'
#' @return data.table
#' @export
#'
#' @examples
tx_add_endRatio1bpUS <- function(DT, minCov = 50){
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("endRatio1bpUS")
    tmp <- (DT$end_3p) / (DT$cov)
    tmp[DT$cov < minCov] <- NA
    DTL <- tibble::add_column(DT, endRatio1bpUS = tmp) %>% txtools::tx_split_DT()
    DT <- lapply(DTL, function(DT){
        DT$endRatio1bpUS <- c(NA, utils::head(DT$endRatio1bpUS, -1))
        DT
    }) %>% txtools::tx_merge_DT()
}

#' Add position unique names
#'
#' Adds a column to DT which pastes the name of the gene with its transcript coordinate
#' creating unique position identifiers, for use in downstream analysis which
#' requires them. The column is added after the mandatory 'txcoor' column.
#' Unique names formation is expected and checked by default.
#'
#' @param DT data.table.
#' @param sep character. Separator between gene and txcoor, by deafult colon sign.
#' @param check_uniq logical. Set to false to override unique position names check.
#'
#' @return data.table
#' @export
#'
#' @author M.A. Garcia-Campos
#' @examples
tx_add_pos <- function(DT, sep = ":", check_uniq = T){
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("pos")
    pos <- paste(DT$gene, DT$txcoor, sep = sep)
    if(check_uniq){
        if(!all(!duplicated(pos))){
            stop("Combinations of gene and txcoor by row in DT are not unique.")
        }
    }
    return(tibble::add_column(DT, pos = pos, .after = "txcoor"))
}

#' Add 1bp-site logical annotation
#'
#' Add a logical variable column in DT, for which the genomic coordinates of a
#' GRanges object are used to differentiate between 1-bp length sites.
#' In which TRUE equals to presence of site in the provided GRanges object.
#' This annotation is useful, for example, when marking already known RNA
#' modification sites.
#'
#' @param DT data.table. Output of txtools data.tables with coverage information
#' as output of tx_coverageDT() and tx_covNucFreqDT()
#' @param GRanges GRanges. Ranges of length 1, to be marked in the data.table
#' @param colName character. Name of the new column to be added.
#' @param nCores integer. Number of cores to use to run function. Multicore
#' capability not available in Windows OS.
#'
#' @return data.table
#' @export
#'
#' @examples
tx_add_siteAnnotation <- function (DT, GRanges, colName, nCores = 1){
    check_mc_windows(nCores)
    check_integerGreaterThanZero_arg(nCores, "nCores")
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent(colName)
    if(!all(GenomicRanges::start(GRanges) == GenomicRanges::end(GRanges))){
        stop("start and ends are not the same in GRanges, only 1-nuc-long sites allowed")
    }
    if(class(GRanges) != "GRanges"){
        stop("GRanges must be of class GRanges")
    }
    DTL <- tx_split_DT(DT)
    tx_merge_DT(parallel::mclapply(mc.cores = nCores, DTL, function(x){
        hlp_add_siteAnnotation(x, GRanges, colName)
    }))
}


#' Adds the presence of a motif
#'
#' @param DT data.table. Output of txtools data.tables with coverage information
#' as output of \code{\link{tx_coverageDT}} \code{\link{tx_covNucFreqDT}}
#' functions.
#' @param motif character. A word which depicts a DNA sequence motif to annotate.
#' The function allows IUPAC ambiguity codes, e.g. R={A|G}.
#' @param nucPositions character or numeric.
#' \itemize{
#' \item "all": Will mark all the positions of the motif
#' \item "center": Will mark the center of the motif
#' \item "i": A number can be passed in which the \strong{i}th position will be
#' marked. For example, for the 'CAC' motif a value of '2' will mark all 'A's
#' surrounded by a 'C'.
#' }
#' @param motifColName character. Name of the new column to be added for
#' annotating motif presence. Automatically is set to be a combination of the
#' input motif and nucPositions arguments.
#' @param mask_N logical. If set to FALSE, 'N' nucleotides are left as is,
#' therefore matching motifs. i.e. A consecutive sequence of NNNNN will match
#' any 5 letter motif as 'DRACH'; generally not desired.
#' @param nCores integer. Number of cores to use to run function. Multicore
#' capability is not available in Windows OS.
#'
#' @return data.table
#' @export
#'
tx_add_motifPresence <- function (DT, motif, nucPositions = "all",
                                  motifColName = "auto", mask_N = TRUE, nCores = 1){
    check_mc_windows(nCores)
    check_integerGreaterThanZero_arg(nCores, "nCores")
    DT <- check_DT(DT)
    if(!("refSeq" %in% names(DT))){
        stop("DT must contain the refSeq column, as added by tx_add_refSeqDT() function.")
    }
    if(mask_N){
        oriSeq <- DT$refSeq
        DT$refSeq[DT$refSeq == "N"] <- "."
    }
    # Process motif
    midMot <- NULL
    if(motifColName == "auto"){
        motifColName <- paste(motif, "motif", paste(nucPositions, collapse = "_"), sep = "_")
    }
    DT <- hlp_removeColumnIfPresent(DT, motifColName)
    if(nucPositions == "center"){
        if((nchar(motif) %% 2) == 0){
            stop("motif has even number of characters, center cannot be determined.")
        }
        midMot <- ceiling(nchar(motif) / 2)
    }else if(is.numeric(nucPositions)){
        if(max(nucPositions) > nchar(motif)){
            stop("nucPositions are bigger than motif length.")
        }
    }else{
        stop("nucPositions must be either: 'all', 'center', or a numeric vector
                 with positions in motif to be marked.")
    }
    oNames <- names(DT)
    DT <- tx_merge_DT(parallel::mclapply(mc.cores = nCores, tx_split_DT(DT), function(x){
        hlp_add_motifPresence(x, motif, nucPositions, midMot)
    })) %>% magrittr::set_names(c(oNames, motifColName))
    if(mask_N){
        DT[, "refSeq"] <- oriSeq
    }
    DT
}

#' Adding rolling mean to txDT
#'
#' @param DT data.table. Output of txtools data.tables with coverage information
#' as output of \code{\link{tx_coverageDT}} \code{\link{tx_covNucFreqDT}}
#' functions.
#' @param colName character. Column to which rolling mean will be applied
#' @param winSize numeric. Size of window
#' @param newColName character. Name to be given to the output. If left empty,
#' will assign a default name, combination of the `colName` and `winSize` arguments.
#' @param fill Either an empty vector (no fill), or a vector (recycled to)
#' length 3 giving left, middle and right fills. By default NA.
#' @param align character. Align windows on the "left", "center" or "right".
#' @param minCov numeric. Minimum coverage required to output ratio. If coverage
#' is less then an NA is output in that position. To output all positions assign 0.
#' @param nCores integer. Number of cores to use to run function. Multicore
#' capability is not available in Windows OS.
#'
#' @return
#' @export
#'
tx_add_rollingMean <- function(DT, colName, winSize, newColName = NULL,
                               fill = NA, align = "center", minCov = 20, nCores = 1){
    if(is.null(newColName)){newColName <- paste("rollMean", colName, winSize, sep = "_" )}
    DT <- hlp_removeColumnIfPresent(DT, newColName)
    if(!is.numeric(DT[[colName]])) stop("Column '", colName, "' must be numeric")
    oNames <- colnames(DT)
    tmp <- parallel::mclapply(mc.cores = nCores, tx_split_DT(DT), function(x){
        RcppRoll::roll_mean(x[[colName]], n = winSize, fill = fill, align = align)
    }) %>% unlist() %>% unname()
    tmp[DT$cov <= minCov] <- NA
    tibble::add_column(DT, tmp) %>% magrittr::set_names(c(oNames, newColName))
}

# Get functions ################################################################

#' Get data from a position and their neighboring positions
#'
#' Get the data from positions marked with TRUE in a logical variable as the
#' center and the values in its vicinity, delimited by a flank size.
#' Useful to aggregate data surrounding specific sites as motif locations.
#'
#' @param DT data.table. Output of txtools data.tables with coverage information
#' as output of \code{\link{tx_coverageDT}} \code{\link{tx_covNucFreqDT}}
#' @param logi_col
#' @param values_col
#' @param upFlank
#' @param doFlank
#' @param addRowNames
#' functions.
#' @return matrix
#' @export
#'
#' @examples
tx_get_flanksFromLogicAnnot <- function(DT, logi_col, values_col, upFlank, doFlank, addRowNames = TRUE){
    if(!is.logical(DT[[logi_col]])){stop("column ", logi_col, " must be of class 'logical'.")}
    if(!values_col %in% colnames(DT)){stop(values_col, " column must be in DT")}
    if(!"gene" %in% colnames(DT)){stop("gene column must be in DT.")}
    Igenes <- as.character(unique(DT[DT[[logi_col]],][["gene"]]))
    DT <- DT[DT$gene %in% Igenes,]
    tmpVal <- split(DT[[values_col]], forcats::fct_drop(DT$gene))
    tmpVar <- split(DT[[logi_col]], forcats::fct_drop(DT$gene))
    spacer <- rep(".", max(upFlank, doFlank))
    spacerVar <- rep(NA, max(upFlank, doFlank))
    fullVal <- c(spacerVar, lapply(tmpVal, function(x) (c(x, spacerVar))) %>% do.call(what = "c"))
    fullVar <- c(spacerVar, lapply(tmpVar, function(x) (c(x, spacerVar))) %>% do.call(what = "c"))
    tmpO <- get_VectorElements(fullVal, which(fullVar) - upFlank, which(fullVar) + doFlank) %>%
        do.call(what = "rbind")
    tickNames <- paste((-upFlank):(doFlank), "bp")
    tickNames[upFlank + 1] <- logi_col
    colnames(tmpO) <- tickNames
    if(addRowNames){
        rownames(tmpO) <- paste(DT$gene[DT[[logi_col]]], DT$txcoor[DT[[logi_col]]], sep = ":")
    }
    return(tmpO)
}

#' Get flanking sequences
#'
#' @param DT
#' @param logi_col
#' @param upFlank
#' @param doFlank
#' @param addNames
#'
#' @return
#' @export
#'
#' @examples
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
        warning("Some sequences reached the end of transcript, a '.' ",
                "was added in place, which may affect downstream results.")
    }
    if(addNames){
        names(tmpSeq) <- paste(DT$gene[DT[[logi_col]]], DT$txcoor[DT[[logi_col]]], sep = ":")
    }
    return(tmpSeq)
}


#' Get length of genes
#'
#' Outputs the length of the transcripts in the DT
#'
#' @param DT data.table. Output of txtools data.tables with coverage information
#' as output of \code{\link{tx_coverageDT}} \code{\link{tx_covNucFreqDT}}
#' functions.
#'
#' @return named numeric. Length of genes as per times they appear in DT
#' @export
#'
#' @examples
tx_get_geneLengths <- function(DT){
    tmpTb <- table(as.character(DT$gene))
    tmpTb %>% as.numeric() %>% magrittr::set_names(names(tmpTb))
}

#' Get metagene at CDS
#'
#' @param txDT
#' @param geneAnnotation
#' @param colVars
#' @param CDS_align
#' @param upFlank
#' @param doFlank
#'
#' @return
#' @export
#'
#' @examples
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
        txDT <- tx_split_DT(txDT) %>% annot_CDSsta_DTL(CDS_start = CDS_start) %>% tx_merge_DT()
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
        txDT <- tx_split_DT(txDT) %>% annot_CDSend_DTL(CDS_end = CDS_end) %>% tx_merge_DT()
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



#' Get transcriptome sequences
#'
#' @param genome
#' @param geneAnnot
#' @param outFile
#' @param nCores
#'
#' @return
#' @export
#'
#' @examples
tx_get_transcriptSeqs <- function(genome, geneAnnot, outFile = NULL, nCores = 1){
    check_GA_genome_chrCompat(geneAnnot, genome)
    CHRS <- as.character(unique(GenomicRanges::seqnames(geneAnnot)))
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


# Statistical tests ############################################################

#' Likelihood Ratio Test
#'
#' @param DTL
#' @param tVar
#' @param sVar
#' @param test_groups
#' @param minTrials
#'
#' @return
#' @export
#'
#' @examples
tx_test_LRTedgeR <- function(DTL, tVar, sVar, test_groups, minTrials = 50){
    DTL <- tx_unifyTxDTL(DTL)
    cMat_cov <- lapply(DTL, function(DT){
        DT[[tVar]]
    }) %>% do.call(what = cbind) %>%
        set_colnames(paste0(names(DTL), "-cov")) %>%
        set_rownames(1:nrow(DTL[[1]]))
    cMat_sta <- lapply(DTL, function(DT){
        DT[[sVar]]
    }) %>% do.call(what = cbind) %>%
        set_colnames(paste0(names(DTL), "-sta")) %>% set_rownames(1:nrow(DTL[[1]]))
    cMat_cov <- cMat_cov - cMat_sta # Difference to make starts + cov = total cov (needed by edgeR)
    # DGE object creation
    Treat <- factor(lapply(test_groups, function(x) rep(x, 2)) %>% unlist)
    smpNames <- paste(lapply(names(DTL), function(x) rep(x, 2)) %>% unlist(),
                      c("-sta", "-cov"), sep = "")
    countTable <- cbind(cMat_sta, cMat_cov)[,smpNames]
    hasNo_NAs <- !rowSums(is.na(countTable)) > 0
    Coverage <- cMat_sta + cMat_cov
    HasCoverage <- rowSums(Coverage >= minTrials) == ncol(Coverage)
    HasBoth <- rowSums(cMat_sta) > 0 & rowSums(cMat_cov) > 0
    cat("Attempting to run in", sum(hasNo_NAs & HasCoverage & HasBoth), "sites\n")
    yall <- edgeR::DGEList(counts = countTable[hasNo_NAs & HasCoverage & HasBoth,], group = Treat)
    yall$genes <- data.frame(DTL[[1]][hasNo_NAs & HasCoverage & HasBoth,c("gene", "txcoor")])
    # Filter by coverage and change in start and coverage (not always being 0)
    RTstoppage <- gl(n = 2, k = 1, length = ncol(yall), labels=c("sta","cov"))
    y <- yall
    # library size
    TotalLibSize <- y$samples$lib.size[RTstoppage=="sta"] + y$samples$lib.size[RTstoppage=="cov"]
    y$samples$lib.size <- rep(TotalLibSize, each=2)
    META <- data.table(names(DTL), group = test_groups)
    # design matrix
    designSL <- model.matrix(~0+group, data = META)
    design <- edgeR::modelMatrixMeth(designSL)
    # Estimate dispersion
    cat("Estimating data dispersion and fitting GNB model\n")
    y1 <- edgeR::estimateDisp(y, design=design, trend = "none", )
    # Fitting Negative Binomial Generalized Linear Models
    fit <- edgeR::glmFit(y1, design)
    #Contrast
    contr <- limma::makeContrasts(onlyIP = (groupWT-groupKO), levels = design)
    #Likelihood ratio test
    lrt <- edgeR::glmLRT(fit, contrast = contr[, "onlyIP"])
    #Results table
    RES_sta <- cbind(lrt$table, y$genes)
    RES_sta$FDR <- p.adjust(RES_sta$PValue, method = "fdr")
    test_mat <- do.call(what = "cbind", lapply(DTL, function(DT) DT[[sVar]] / DT[[tVar]]))
    rownames(test_mat) <- paste(DTL[[1]][["gene"]], DTL[[1]][["txcoor"]], sep = ":")
    test_mat <- test_mat[paste(y$genes$gene, y$genes$txcoor, sep = ":"),] %>%
        rowMeansColG(test_groups)
    test_mat <- test_mat[,1] - test_mat[,2]
    RES_sta <- data.table(RES_sta, RD = test_mat)
    RES_sta <- RES_sta[order(RES_sta$FDR, decreasing = FALSE),]
    if("refSeq" %in% colnames(DTL[[1]])){
        RES_sta <- dplyr::right_join(RES_sta,
                                     select(DTL[[1]],
                                            c("chr", "gencoor", "strand",
                                              "gene", "txcoor", "refSeq")),
                                     by = c("gene", "txcoor"))
    }else{
        RES_sta <- dplyr::right_join(RES_sta,
                                     select(DTL[[1]],
                                            c("chr", "gencoor", "strand",
                                              "gene", "txcoor")))
    }
    return(RES_sta)
}


#' t test
#'
#' @param DTL
#' @param cont_var
#' @param test_groups
#' @param test_na.rm
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
tx_test_ttest <- function(DTL, cont_var, test_groups, test_na.rm = FALSE, ...){
    sapply(DTL, function(x) check_DThasCol(x, cont_var))
    #check if needed to unify
    if(!check_sameGenesInDTL(DTL)){
        cat("Using the intersection of genes on DTL")
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

# Other accessory functions ####################################################

#' Centered numeric sequence
#'
#' Creates a numerical sequence which is centered in the first argument and
#' is of length twice the second argument plus one.
#'
#' @param position numeric
#' @param windowLength numeric
#'
#' @return numeric
#' @export
#'
#' @examples
#' # Numeric interval centered in 10, with 5 closest negative
#' # and positive integers
#' window_around(10, 5)
window_around <- function(position, windowLength){
    (position - windowLength):(position + windowLength)
}

#' Total counts of reads per gene model
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via the tx_reads() function.
#'
#' @return integer
#' @export
#'
#' @author M.A. Garcia-Campos
#'
#' @examples
tx_counts <- function(x){
    tmp <- suppressWarnings(unlist(x)) %>% GenomeInfoDb::seqnames() %>% table()
    magrittr::set_names(as.integer(tmp), names(tmp))
}

# Built-in data loading functions ##############################################

#' D. melanogaster gene annotation subset path
#'
#' Returns the local path to the in-package BED file containing gene models
#' for highly expressed genes in the 'Pasilla' experiment.
#'
#' @return character. Path to dm3 gene annotation file
#' @export
#'
#' @examples
tx_dm3_geneAnnot <- function(){
    system.file("extdata", "toyGeneAnnot_Dmelan_chr4.bed", package = "txtools")
}
