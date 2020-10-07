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
#' # Loading in-package BAM file
#' bamFile <- system.file("extdata", "example_hg19.bam", package = "txtools")
#' hg19_bam <- tx_load_bam(bamFile, pairedEnd = TRUE, loadSeq = TRUE, verbose = TRUE)
#' summary(hg19_bam)
tx_load_bam <- function(file,
                        pairedEnd,
                        yieldSize = 100000,
                        scanFlag = "default",
                        loadSeq = FALSE,
                        recoverDumpedAligns = FALSE,
                        verbose = TRUE){
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
    if(length(readCycles) > 200){
        warning("If taking too much time, try to subset your BAM file. \n")
    }
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
    bamData <- list(GAligns = reads,
                    dumpedAmbigPairs = GenomicAlignments::getDumpedAlignments())
    if(verbose){cat(" \n")}
    if(verbose){
        cat(length(bamData$GAligns), "Reads succesfully loaded \n")
        cat("Dumped ambioguous reads:", length(bamData$dumpedAmbigPairs), "\n")
    }
    if(recoverDumpedAligns == F){
        if(length(bamData$dumpedAmbigPairs) > 0){
            warning(length(bamData$dumpedAmbigPairs), " dumped ambiguous reads. ",
                    "Use 'GenomicAlignments::getDumpedAlignments()' to retrieve ",
                    "them from the dump environment, or set the ",
                    "'recoverDumpedAligns' argument to TRUE")
        }
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

#' Filter ranges by a maximum width
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via the tx_reads() or tx_reads_mc()
#' functions.
#' @param thr integer. Threshold for maximum width size allowed on output.
#' @param nCores integer. Number of cores to use to run function. Multicore
#' capability not available in Windows OS.
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

#' Transcriptomic reads convertion
#'
#' Assign aligned reads to their respective gene models and convert their
#' positions into a transcriptomic coordinate system. It also stitches together
#' paired-end aligned reads into a single character 'word' in which dots '.'
#' separate Read1 and Read2 by their corresponding insert.
#'
#' @param reads GAlignments or GAlignmentPairs. Genomic alignments to be processed
#' @param geneAnnot GenomicRanges. Gene annotation loaded via the tx_load_bed()
#' @param minReads integer. Minimum number of reads required to overlap a gene
#' @param withSeq logical. Set to TRUE if sequence should be preserved; 'reads'
#' object should contain sequences.
#' @param verbose logical. Set to FALSE to show less information.
#' @param nCores integer. Number of cores to use to run function. Multicore
#' capability not available in Windows OS.
#'
#' @return GRanges
#' space. Alignments are located in genes, instead of chromosomes.
#' @export
tx_reads <- function(reads, geneAnnot, minReads = 50, withSeq = F, verbose = T,
                     nCores = 1){
    # Checks
    check_mc_windows(nCores)
    check_integerGreaterThanZero_arg(minReads, "minReads")
    check_GA_reads_compatibility(reads, geneAnnot)
    check_BAM_has_seq(reads)
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
            cat(length(unique(unlist(split_i))), "reads overlap",
                length(geneAnnot), "gene models \n")
            cat("Filtering reads by gene model... \n")
            if(withSeq){
                cat("Processing sequences. This may take several minutes... \n")
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
    return(OUT)
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
#' @author M.A. Garcia-Campos
#' @aliases tx_coverageDT
#'
#' @examples
tx_makeDT_coverage <- function(x, geneAnnot, genome = NULL, fullDT = FALSE,
                               nCores = 1){
    check_mc_windows(nCores)
    check_integerGreaterThanZero_arg(nCores, "nCores")
    OUT <- hlp_cbind2Tabs(hlp_genCoorTab_mc(x, geneAnnot, nCores),
                       hlp_coverageTab_mc(x, nCores)) %>% tx_merge_DT()
    if(fullDT){
        OUT <- tx_complete_DT(DT = OUT, geneAnnot = geneAnnot, genome = genome, nCores = nCores)
    }
    if(!is.null(genome)){
        OUT <- tx_add_refSeqDT(OUT, genome = genome, geneAnnot = geneAnnot, nCores = nCores)
    }else{return(OUT)}
}

#' Summarized Nucleotide Frequency data.table
#'
#' This function constructs a list of data.tables that contains nucleotide frequency
#' metrics per nucleotide by transcript:
#'\itemize{
#'\item cov = Insert coverage
#'\item start_5p = read-start counts
#'\item end_3p = read-end counts
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
#' @author M.A. Garcia-Campos
#' @aliases tx_nucFreqDT
#'
#' @examples
tx_makeDT_nucFreq <- function(x, geneAnnot, genome = NULL,
                              simplify_IUPAC = "splitForceInt", fullDT = FALSE,
                              nCores = 1){
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
#' alignments data by gene. Constructed via the tx_reads() or tx_reads_mc()
#' functions.
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
#' @author M.A. Garcia-Campos
#' @aliases tx_covNucFreqDT
#'
#' @examples
tx_makeDT_covNucFreq <- function(x, geneAnnot, genome = NULL,
                                 simplify_IUPAC = "splitForceInt", fullDT = FALSE,
                                 nCores = 1){
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
    tmp <- split(DT, by = "gene", drop = dropEmpty)
    lapply(tmp, function(y) {
        y[order(y$txcoor), ]
    })
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
    tmpCoorTabs <- hlpr_genCoorTabGenes(missGenes, geneAnnot, genome, nCores) %>%
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
    tx_merge_DT(completeDTL)
}

# data.table functions #########################################################

#' Add reference sequence to data.table
#'
#' @param DT data.table. A summarized data.table object. See tx_coverageDT(),
#' tx_nucFreq() and tx_covNucFreqDT() functions.
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
    check_mc_windows(nCores)
    check_integerGreaterThanZero_arg(nCores, "nCores")
    DT <- check_DT(DT)
    check_GA_genome_chrCompat(geneAnnot = geneAnnot, genome = genome)
    DTL <- tx_split_DT(DT)
    tmp <- parallel::mclapply(mc.cores = nCores, DTL, function(DT){
        hlp_add_refSeqDT(DT, fastaGenome = genome, geneAnnot = geneAnnot)
    })
    OUT <- tx_merge_DT(tmp)
    return(OUT)
}

# Add column of sum of nucleotide frequency different to the reference sequence
# Counting deletions, not counting 'N's and inserts into calculation.
tx_add_diffNucToRef <- function(DT){
    DT <- check_DT(DT)
    selNucs <- setdiff(intersect(txtools::IUPAC_code_2nucs, names(DT)), c(".", "N"))
    tmp <- DT[, selNucs, with = FALSE]
    OUT <- rep(NA, nrow(DT))
    nucsInRef <- unique(DT$refSeq)
    for(i in nucsInRef){
        selNucs <- setdiff(c("A", "C", "G", "T", "-"), i)
        OUT[which(DT$refSeq == i)] <- rowSums(tmp[which(DT$refSeq == i), selNucs, with = F])
    }
    tibble::add_column(DT, diffToRef = OUT)
}

# Add sum of nucleotide frequencies for meaningful nucleotides (i.e. not counting Ns and inserts)
tx_add_nucTotal <- function(DT){
    DT <- check_DT(DT)
    selNucs <- setdiff(intersect(txtools::IUPAC_code_2nucs, names(DT)), c(".", "N"))
    out <- rowSums(DT[, selNucs, with = FALSE])
    tibble::add_column(DT, nucTotal = out)
}

# Add column of Different Nucleotide to reference ratio, diffToRef and nucTotal columns are required
tx_add_diffNucToRefRatio <- function(DT, addDiffandTotalCols = FALSE){
    DT <- check_DT(DT)
    tmpDT <- tx_add_diffNucToRef(DT) %>% tx_add_nucTotal()
    tmp <- round(tmpDT$diffToRef / tmpDT$nucTotal, 6)
    if(addDiffandTotalCols){DT <- tmpDT}
    tibble::add_column(DT, diffToRefRatio = tmp)
}

#' Add startRatio to data.table
#'
#' Add read-starts ratio over coverage.
#'
#' @param DT data.table. Output of txtools data.tables with coverage information
#' as output of \code{\link{tx_coverageDT}} \code{\link{tx_covNucFreqDT}}
#' functions
#'
#' @return data.table
#' @export
#'
tx_add_startRatio <- function(DT){
    DT <- check_DT(DT)
    tibble::add_column(DT, startRatio = DT$start_5p / DT$cov, .after = "start_5p")
}

#' Add endRatio to data.table
#'
#' Add read-ends ratio over coverage.
#'
#' @param DT data.table. Output of txtools data.tables with coverage information
#' as output of \code{\link{tx_coverageDT}} \code{\link{tx_covNucFreqDT}}
#' functions
#'
#' @return data.table
#' @export
#'
tx_add_endRatio <- function(DT){
    DT <- check_DT(DT)
    tibble::add_column(DT, endRatio = DT$end_3p / DT$cov, .after = "end_3p")
}

#' Add position names column to DT
#'
#' Adds a column which pastes the name of the gene with its transcript coordinate
#' creating a unique position identifier, for use in downstream analysis which
#' requires it. The column is added after the 'txcoor' column. Unique names are
#' checked as well.
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
    DT <- check_DT(DT)
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
#' GRanges object are used to mark differentiate between sites in the GRanges
#' object and all others. TRUE = site in GRanges object. This is useful for
#' example when marking already known RNA modification sites.
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
    DT <- check_DT(DT)
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

# Other accesory functions #####################################################

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
