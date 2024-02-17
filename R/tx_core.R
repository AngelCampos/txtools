#' txtools: A package facilitating analysis of RNA modifications, structures, and interactions
#'
#' txtools enables the processing, analysis, and visualization of RNA-seq data
#' at the nucleotide-level resolution, seamlessly integrating alignments to
#' the genome with transcriptomic representation. txtools’ main inputs are
#' BAM files and a transcriptome annotation, and the main output is a table,
#' capturing mismatches,  deletions, and the number of reads beginning and
#' ending at each nucleotide in the transcriptomic space. txtools further
#' facilitates downstream visualization and analyses.
#'
#' Most of txtools' functions start with the prefix tx_ and are grouped by
#' families:
#'
#' @section tx_load_*():
#' Load initial data as genomes (FASTA), gene annotations (BED), and mapped
#' reads (BAM).
#' @section tx_add_*():
#' Add a new variable to the txDT, generally by computing a ratio or frequency.
#' Their output is the new txDT. e.g. \code{\link{tx_add_startRatio}}(), which
#' adds the start to coverage ratio; \code{\link{tx_add_motifPresence}}(),
#' which adds the location of RNA sequence motifs across the transcriptome
#' @section tx_get_*():
#' Extract information from a txDT and generate an object that is NOT a txDT.
#' e.g. \code{\link{tx_get_metageneRegions}}() which outputs a metagene matrix
#' with each row representing a gene and each column a bin in one of the
#' codifying gene regions
#' @section tx_plot_*():
#' Plotting functions. e.g. \code{\link{tx_plot_nucFreq}}() and
#' \code{\link{tx_plot_staEndCov}}(), which plot the counts of data of
#' nucleotide frequency, and read-starts/ends and coverage respectively.
#' @section tx_test_*():
#' Use of the txDT objects from experimental data to do statistical tests of
#' metrics between groups of samples. e.g. \code{\link{tx_test_ttest}}() which
#' performs t-tests using a list of txDTs and a vector of the groups.
#'
#' @docType package
#' @name txtools
#' @aliases txtools-package
#' @keywords DataImport DataRepresentation Coverage Epitranscriptomics RNASeq
#' @keywords Transcription SNP
#' @importFrom rlang .data
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
#' @param strandMode numeric. \itemize{
#' \item 1 (default): Strand of the pair is that of its **first** alignment: Directional Illumina (Ligation), Standard SOLiD. (Single-end No change in strand)
#' \item 2: strand of the pair is strand of its **last** alignment: dUTP, NSR, NNSR, Illumina stranded TruSeq PE protocol. (Single-end: Change to inverse strand)
#' \item 0: strand of the pair is set to '*' (unspecified). This mode is no longer supported by \code{\link{tx_reads}}() .} More info: \code{\link[GenomicAlignments]{GAlignmentPairs-class}}.
#' @param verbose logical. Set to FALSE to show less information.
#'
#' @return GRanges
#' @export
#'
#' @examples
#' # Loading in-package BAM file
#' bamFile <- system.file("extdata", "example_hg19.bam", package = "txtools")
#' hg19_bam <- tx_load_bam(bamFile, pairedEnd = TRUE, loadSeq = TRUE, verbose = TRUE)
#' summary(hg19_bam)
tx_load_bam <- function(file, pairedEnd, yieldSize = 100000,
                        scanFlag = "default", loadSeq = FALSE, strandMode = 1,
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
    #  Opens and reads BAM files in chunks of length = yieldSize
    open(BAMFILE)
    reads <- lapply(pbJumps, function(i){
        if(verbose){utils::setTxtProgressBar(pb, i)}
        if(pairedEnd){
            suppressWarnings(
                GenomicAlignments::readGAlignmentPairs(BAMFILE,
                                                       use.names = TRUE,
                                                       param = readParams,
                                                       strandMode = strandMode))
        }else if(pairedEnd == FALSE){
            suppressWarnings(
                GenomicAlignments::readGAlignments(BAMFILE,
                                                   use.names = TRUE,
                                                   param = readParams))
        }
    }) %>% do.call(what = "c")
    close(BAMFILE)
    if(verbose){close(pb)}
    if(pairedEnd == FALSE & strandMode == 2){
        reads <- BiocGenerics::invertStrand(reads)
    }else if(pairedEnd == FALSE & strandMode == 3){
        GenomicAlignments::strand(reads) <- "*"
    }
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
tx_load_bed <- function(bedfile){
    tmp <- plyranges::read_bed(bedfile)
    if(length(S4Vectors::mcols(tmp)) == 2){
        tmp$itemRGgb <- NA
        tmp$thick <- tmp@ranges
        tmp$blocks <- IRanges::IRangesList(S4Vectors::splitAsList(
            IRanges::IRanges(start = 1, end = IRanges::end(tmp) - IRanges::start(tmp) + 1),
            factor(seq(length(tmp)))), compress = TRUE) %>%
            magrittr::set_names(NULL)
    }
    # Check no duplicated gene names
    dupN <- duplicated(tmp$name)
    if(sum(dupN) > 0){
        stop("Duplicated genes found in gene annotation.\n",
             paste(tmp$name[dupN], collapse = " "))
    }
    return(tmp)
}

#' Load genome
#'
#' Load genome as DNAStrinSet from FASTA file.
#'
#' Chromosome names will be limited to the first word separated by a space
#'
#' @param fastaFile path to FASTA file with genomic sequences
#'
#' @return DNAStringSet
#' @export
tx_load_genome <- function(fastaFile){
    out <- Biostrings::readDNAStringSet(fastaFile)
    names(out) <- names(out) %>%
        stringr::str_split(pattern = " ") %>%
        lapply(function(x) x[1]) %>% unlist()
    out
}

#' Loading RDS files into data.tables
#'
#' @param file File name or path to .rds file to be loaded
#'
#' @return data.table
#' @export
#'
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
#' To retrieve unassigned alignments use the function tx_get_unassignedAlignments()
#'
#' @param reads GAlignments or GAlignmentPairs. Genomic alignments to be processed
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}()
#' @param minReads integer. Minimum number of alignments required to overlap a gene
#' @param withSeq logical. Set to TRUE if sequence should be preserved; 'reads'
#' object should contain sequences.
#' @param verbose logical. Set to FALSE to show less information.
#' @param nCores integer. Number of cores to run the function with. Multi-core
#' capability not available in Windows OS.
#' @param ignore.strand logical. Set to TRUE to allow alignments in a gene
#' ignoring the strand of the alignment. False by default.
#'
#' @aliases tx_reads_mc
#' @aliases tx_flushUnassigned
#'
#' @return GRanges
#' @export
tx_reads <- function(reads, geneAnnot, minReads = 50, withSeq = FALSE,
                     verbose = TRUE, ignore.strand = FALSE, nCores = 1){
    # Checks
    tx_flushUnassigned()
    check_mc_windows(nCores)
    check_integerGreaterThanZero_arg(minReads, "minReads")
    check_GA_reads_compatibility(reads, geneAnnot)
    if(withSeq){check_BAM_has_seq(reads)}
    if(!class(reads) %in% c("GAlignmentPairs", "GAlignments")){
        stop("reads argument should be of class GAlignmentPairs or GAlignments. \n")
    }
    reads <- rmAmbigStrandAligns(reads)
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
                      GAlignmentPairs = hlpr_splitReadsByGenes(
                          reads = reads, bedR = geneAnnot, overlapType = overlapType,
                          minReads = minReads, ignore.strand = ignore.strand),
                      GAlignments = hlpr_splitReadsByGenes_singleEnd(
                          reads = reads, bedR = geneAnnot, overlapType = overlapType,
                          minReads = minReads, ignore.strand = ignore.strand))
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
    }else{stop("There where no alignments, or not enough alignments overlapping the gene models. \n")}
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
                             minReads = minReads,
                             ignore.strand = ignore.strand)
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
    OUT <- OUT[lapply(OUT, length) %>% unlist %>% magrittr::is_weakly_greater_than(minReads)] %>%
        GenomicRanges::GRangesList()
    if(verbose){
        cat("Output contains:", lapply(OUT, names) %>% unlist %>% unique %>% length,
            "unique alignments in", length(OUT), "gene models \n")
    }
    # Dump not assigned alignments into dump environment
    if(length(setdiff(names(reads), unique(unlist(lapply(OUT, names))))) > 1L){
        warning("Some alignments were not assigned to any gene, you can ",
                "retrieve them using the tx_get_unassignedAlignments() function.")
        hlp_dump_notAssigned(reads, OUT)
    }
    return(OUT)
}

#' Retrieve dumped alignments
#'
#' Retrieve the read mappings that were not assigned to any gene by
#' \code{\link{tx_reads}}
#'
#' @return GAlignments
#' @export
#'
#' @aliases tx_getUnassignedAlignments
tx_get_unassignedAlignments <- function(){
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
tx_extend_UTR <- function(GR, ext_5p = 0, ext_3p = 0){
    if(is.null(GR$blocks)){
        stop("GR does not containg 'blocks' meta column.")
    }
    if(class(GR$blocks) == "IRanges"){
        GR$blocks <- IRanges::IRangesList(start = split(rep(1, length(GR)), 1:length(GR)),
                                          end = split(rep(1, length(GR)), 1:length(GR)))
    }
    GR_out <- plyranges::stretch(plyranges::anchor_5p(GR), ext_3p)
    GR_out <- plyranges::stretch(plyranges::anchor_3p(GR_out), ext_5p)
    GR_out$blocks <- stretchBlocks_3p(GR_out$blocks, ext_3p, GenomicRanges::strand(GR))
    GR_out$blocks <- stretchBlocks_5p(GR_out$blocks, ext_5p, GenomicRanges::strand(GR))
    return(GR_out)
}

#' Cut gene annotation by txDT's genes
#'
#' @param txDT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}()
#'
#' @return data.table
#' @export
tx_cut_geneAnnotBytxDT <- function(txDT, geneAnnot){
    geneAnnot[match(unique(txDT$gene), geneAnnot$name)]
}

# Manipulate GRangesList #######################################################

#' Filter ranges by a maximum width
#'
#' @param x CompressedGRangesList. Genomic Ranges list containing genomic
#' alignments data by gene. Constructed via tx_reads().
#' @param thr integer. Threshold for maximum width size allowed on output.
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability is not available in Windows OS.
#'
#' @return CompressedGRangesList
#' @export
#'
#' @author M.A. Garcia-Campos
#' @aliases tx_filter_max_width
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
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability not available in Windows OS.
#'
#' @return CompressedGRangesList.
#' @export
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
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}()
#' @param genome list. The full reference genome sequences, as loaded by
#' \code{\link{tx_load_genome}}() or prepackaged by BSgenome, see ?BSgenome::available.genomes
#' @param fullDT logical. Set to TRUE if it is desired to output a data.table
#' with all genes and in the same order as 'geneAnnot' object.
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability not available in Windows OS.
#' @return data.table
#' @export
#'
#' @family makeDT functions
#' @author M.A. Garcia-Campos
#' @aliases tx_coverageDT
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
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}()
#' @param genome list. The full reference genome sequences, as loaded by
#' \code{\link{tx_load_genome}}() or prepackaged by BSgenome, see ?BSgenome::available.genomes
#' @param simplify_IUPAC string. Available options are :
#' \itemize{
#' \item "not": Will output the complete nucleotide frequency table including
#' ambiguous reads using the IUPAC ambiguity code.
#' See: \code{\link[Biostrings]{IUPAC_CODE_MAP}}
#' \item "splitForceInt" (Default): Will force an integers split in which ambiguous codes
#' will be split and assigned half the frequency into their respective nucleotides,
#' if the frequency is an odd number the uneven count will be assigned as "N".
#' \item "splitHalf": Ambiguous nucleotide frequencies will be split in half to
#' their corresponding nucleotides, in cases where frequency is odd creating
#' non-integer frequencies.
#' }
#' @param fullDT logical. Set to TRUE if it is desired to output a data.table
#' with all genes and in the same order as 'geneAnnot' object.
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability not available in Windows OS.
#' @return data.table
#' @export
#'
#' @family makeDT functions
#' @author M.A. Garcia-Campos
#' @aliases tx_nucFreqDT
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
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}()
#' function.
#' @param genome list. The full reference genome sequences, as loaded by
#' \code{\link{tx_load_genome}}() or prepackaged by BSgenome, see ?BSgenome::available.genomes
#' @param simplify_IUPAC string. Available options are :
#' \itemize{
#' \item "not": Will output the complete nucleotide frequency table including
#' ambiguous reads using the IUPAC ambiguity code.
#' See: \code{\link[Biostrings]{IUPAC_CODE_MAP}}
#' \item "splitForceInt" (Default): Will force an integers split in which ambiguous codes
#' will be split and assigned half the frequency into their respective nucleotides,
#' if the frequency is an odd number the uneven count will be assigned as "N".
#' \item "splitHalf": Ambiguous nucleotide frequencies will be split in half to
#' their corresponding nucleotides, in cases where frequency is odd creating
#' non-integer frequencies.
#' }
#' @param fullDT logical. Set to TRUE if it is desired to output a data.table
#' with all genes and in the same order as 'geneAnnot' object.
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability not available in Windows OS.
#' @return data.table
#' @export
#'
#' @family makeDT functions
#' @author M.A. Garcia-Campos
#' @aliases tx_covNucFreqDT
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
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param dropEmpty logical. Drops empty list elements, which occurs when data
#' of genes have been entirely removed, but kept listed in the x$gene factor levels.
#'
#' @return list
#' @export
tx_split_DT <- function(DT, dropEmpty = TRUE){
    split(DT, f = DT$gene, drop = dropEmpty)
}

#' Cutting 5' and 3' ends of data.table using txcoors
#'
#' This function removes the heading (5'UTR) and trailing (3' UTR) portions of
#' data.tables.
#'
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param cut_5p integer Basepairs to me removed from the start of DT for each gene
#' @param cut_3p integer Basepairs to me removed from the end of DT for each gene
#'
#' @return data.table
#' @export
tx_cutEnds_DT <- function(DT, cut_5p = 0, cut_3p = 0){
    DTL <- tx_split_DT(check_DT(DT))
    lapply(DTL, function(x){
        hlp_remove_UTR(x, cut_5p, cut_3p)
    }) %>% tx_merge_DT()
}

#' Order txDT
#'
#' Orders a txDT by gene name and transcriptomic coordinate (txcoor).
#' Accordingly reorders chromosome, strand, and gene factors' levels
#'
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#'
#' @return data.table
#' @export
tx_orderDT <- function(DT){
    DT$chr <- factor(as.character(DT$chr))
    DT$strand <- factor(as.character(DT$strand))
    DT$gene <- factor(as.character(DT$gene))
    DT[order(as.character(DT$gene), DT$txcoor), ]
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
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().
#' @param genome list. The full reference genome sequences, as loaded by
#' \code{\link{tx_load_genome}}() or prepackaged by BSgenome, see ?BSgenome::available.genomes
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability not available in Windows OS.
#'
#' @return data.table
#' @export
#'
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
                                        nCores = nCores) %>% tx_merge_DT()
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
                             magrittr::set_names(missCols))
    tx_orderDT(rbind(DT, tmpCoorTabs))
}


#' Shift column in txDT
#'
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param colToShift character. Name of column to be shifted
#' @param direction character. Direction of shift, either 'downstream' or 'upstream'.
#' @param bp integer. Number of nucleotides to shift data to
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability not available in Windows OS.
#' @return data.table
#' @export
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
#' @param txDTL list. A list of txDTs
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().
#' @param genome list. The full reference genome sequences, as loaded by
#' \code{\link{tx_load_genome}}() or prepackaged by BSgenome, see ?BSgenome::available.genomes
#' @param type character. Set function to apply to genes in the list of txDTs, either
#' 'union' or 'intersection'.
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability not available in Windows OS.
#' @return list
#' @export
tx_unifyTxDTL <- function(txDTL, geneAnnot = NULL, genome = NULL, type = "intersection", nCores = 1){
    if(!all(Reduce(intersect, lapply(txDTL, colnames)) == Reduce(union, lapply(txDTL, colnames)))){
        stop("Column names of the elements of txDTL must be the same")
    }
    if(type == "intersection"){
        selGenes <- Reduce(intersect, lapply(txDTL, function(x) unique(x$gene)))
        parallel::mclapply(mc.cores = nCores, txDTL, function(x){
            x <- x[x$gene %in% selGenes,]
            tx_orderDT(x)
        })
    }else if(type == "union"){
        if(is.null(geneAnnot)){stop("geneAnnot must be provided, as loaded with tx_load_bed()")}
        if(is.null(genome)){stop("geneAnnot must be provided, as loaded with tx_load_genome()")}
        selGenes <- Reduce(union, lapply(txDTL, function(x) unique(x$gene)))
        selGA <- geneAnnot[geneAnnot$name %in% selGenes]
        parallel::mclapply(mc.cores = nCores, txDTL, function(x){
            x <- x[x$gene %in% selGenes, ]
            tx_complete_DT(DT = x, geneAnnot = selGA, nCores = nCores, genome = genome)
        })
    }else{stop("Argument 'type' must be either 'intersection' or 'union'")}
}

#' Sample txDT by genes
#'
#' @param txDT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param size number of genes to be sampled
#'
#' @return data.table
#' @export
tx_sampleByGenes <- function(txDT, size){
    GENES <- as.character(unique(txDT$gene))
    if(size > length(GENES)){stop("size cannot be greater than number of genes in txDT")}
    smpGENES <- sample(GENES, size)
    txDT[txDT$gene %in% smpGENES,]
}

# data.table functions #########################################################

#' Add reference sequence to a data.table
#'
#' Adds a column of characters representing the corresponsing nucleotide at
#' the reference genome position, considering the strand of the transcript.
#' The input DT must be the output, or resemble it, of the
#'
#'
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param genome list. The full reference genome sequences, as loaded by
#' \code{\link{tx_load_genome}}() or prepackaged by BSgenome, see ?BSgenome::available.genomes
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability not available in Windows OS.
#'
#' @return data.table
#' @export
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
#' @param DT data.table. A table as output by the
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#'
#' @return data.table
#' @export
tx_add_misincCount <- function(DT){
    check_refSeq(DT)
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("misincCount")
    selNucs <- setdiff(intersect(txtools::IUPAC_code_2nucs, names(DT)), c(".", "N"))
    tmp <- DT[, selNucs, with = FALSE]
    OUT <- rep(NA, nrow(DT))
    nucsInRef <- unique(DT$refSeq)
    for(i in nucsInRef){
        selNucs <- setdiff(c("A", "C", "G", "T", "-"), i)
        OUT[which(DT$refSeq == i)] <- rowSums(tmp[which(DT$refSeq == i), selNucs, with = FALSE])
    }
    tibble::add_column(DT, misincCount = OUT)
}


#' Add total number of nucleotide reads
#'
#' Add a column to DT of the sum of nucleotide frequencies for meaningful
#' nucleotide reads not counting undetermined 'N' and inserts '.'.
#'
#' @param DT data.table. A table as output by the
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#'
#' @return data.table
#' @export
#'
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
#' @param DT data.table. A table as output by the
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param addCounts Set to TRUE to add counts of total nucleotides
#' read (nucTotal) and different to reference nucleotide (misincCount) columns.
#' @param minNucReads Minimum number of nucleotides read needed to calculate
#' misincorporation rate
#'
#' @return data.table
#' @export
#'
#' @seealso \code{\link{tx_add_diffNucToRef}} and \code{\link{tx_add_nucTotal}}
tx_add_misincRate <- function(DT, minNucReads = 20, addCounts = FALSE){
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("misincRate")
    tmpDT <- tx_add_misincCount(DT) %>% tx_add_nucTotal()
    tmp <- round(tmpDT$misincCount / tmpDT$nucTotal, 6)
    tmp[tmpDT$nucTotal < minNucReads] <- NA
    if(addCounts){DT <- tmpDT}
    tibble::add_column(DT, misincRate = tmp)
}
#' Nucleotide specific misincorporation rate
#'
#' Calculates the misincorporation rate from a \bold{ref}erence nucleotide to an
#' specific \bold{mis}incorporated nucleotide. For example, when looking for
#' positions in cytidines that were read as thymine in RNA after bisulphite treatment.
#'
#' @param DT data.table. A table as output by the
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param refNuc Reference nucleotide See \code{\link{tx_add_refSeqDT}}()
#' @param misNuc Misincorporated nucleotide
#' @param minNucReads Minimum number of nucleotides read needed to calculate
#' misincorporation rate
#'
#' @return data.table
#' @export
#'
#' @aliases tx_add_misincorpRateNucSpec
tx_add_misincRateNucSpec <- function(DT, refNuc, misNuc, minNucReads = 20){
    DT <- data.table::data.table(DT)
    newColName <- paste0("MR_", refNuc, "to", misNuc)
    if(newColName %in% colnames(DT)){
        DT <- DT[,!grepl(newColName, colnames(DT))]
    }
    tmp <- round(DT[[misNuc]] / (DT[[misNuc]] + DT[[refNuc]]), 6)
    tmpNucReads <- DT[[misNuc]] + DT[[refNuc]]
    tmp[tmpNucReads < minNucReads] <- NA
    tmp[!(DT$refSeq %in% refNuc)] <- NA
    DT[, newColName] <- tmp
    DT
}

#' Add starts to coverage ratio
#'
#' Add a column to DT of the read-starts to coverage ratio.
#'
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
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
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param minCov numeric. Minimum coverage required to output ratio. If coverage
#' is less then an NA is output in that position.
#'
#' @return data.table
#' @export
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
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param minCov numeric. Minimum coverage required to output ratio. If coverage
#' is less then an NA is output in that position.
#'
#' @return data.table
#' @export
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
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
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
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}()
#' or \code{\link{tx_makeDT_covNucFreq}}() functions.
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
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}()
#' or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param minCov numeric. Minimum coverage required to output ratio. If coverage
#' is less then an NA is output in that position. To output all positions assign 0.
#'
#' @return data.table
#' @export
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
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param sep character. Separator between gene and txcoor, by deafult colon sign.
#' @param check_uniq logical. Set to false to override unique position names check.
#'
#' @return data.table
#' @export
#'
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
#' GRanges object are used to annotate 1 bp sites.
#' In which TRUE equals to presence of the site in the provided GRanges object.
#' This annotation is useful, for example, when marking already known RNA
#' modification sites.
#'
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param GRanges GRanges. Ranges of length 1, to be marked in the data.table
#' @param colName character. Name of the new column to be added.
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability not available in Windows OS.
#'
#' @return data.table
#' @export
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

#' Add motif presence
#'
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
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
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability is not available in Windows OS.
#'
#' @return data.table
#' @export
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

#' Add rolling mean
#'
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param colName character. Column to which rolling mean will be applied
#' @param winSize numeric. Size of window
#' @param newColName character. Name to be given to the output. If left empty,
#' will assign a default name, combination of the `colName` and `winSize` arguments.
#' @param fill Either an empty vector (no fill), or a vector (recycled to)
#' length 3 giving left, middle and right fills. By default NA.
#' @param align character. Align windows on the "left", "center" or "right".
#' @param minCov numeric. Minimum coverage required to output ratio. If coverage
#' is less then an NA is output in that position. To output all positions assign 0.
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability is not available in Windows OS.
#'
#' @return data.table
#' @export
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

#' Add SpliceSites
#'
#' @param txDT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().
#'
#' @return data.table
#' @export
tx_add_spliceSitesLogical <- function(txDT, geneAnnot){
    GENES <- as.character(unique(txDT$gene))
    GA_GR <- exonGRanges(geneAnnot[geneAnnot$name %in% GENES])
    tx_merge_DT(hlp_splLog_v(iGene = GENES, GA_GR, tx_split_DT(txDT)))
}

#' Add gene regions
#'
#' @param txDT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability is not available in Windows OS.
#'
#' @return data.table
#' @export
tx_add_geneRegion <- function(txDT, geneAnnot, nCores = 1){
    check_GA_txDT_compat(txDT, geneAnnot)
    txDT <- hlp_removeColumnIfPresent(txDT, "geneRegion")
    invisible(sapply(c("gencoor", "strand", "gene"), function(x) check_DThasCol(txDT, x)))
    geneAnnot <- geneAnnot[geneAnnot$name %in% unique(txDT$gene)]
    NCG <- geneAnnot[GenomicRanges::width(geneAnnot$thick) == 0]
    CG <- geneAnnot[GenomicRanges::width(geneAnnot$thick) > 0]
    if (length(NCG) > 0) {
        warning(length(NCG), " non-coding genes where marked as 'non-coding'.")
    }
    if (sum(GenomicRanges::strand(CG) == "*") > 0) {
        stop("Genes with no set strand (*) are not allowed in geneAnnot.")
    }
    pos_CG <- as.factor(GenomicRanges::strand(CG)) == "+"
    neg_CG <- as.factor(GenomicRanges::strand(CG)) == "-"
    CDS_tab <- rbind(data.frame(gene = CG[pos_CG]$name,
                                start = IRanges::start(CG[pos_CG]$thick),
                                end = IRanges::end(CG[pos_CG]$thick),
                                strand = "+"),
                     data.frame(gene = CG[neg_CG]$name,
                                start = IRanges::end(CG[neg_CG]$thick),
                                end = IRanges::start(CG[neg_CG]$thick),
                                strand = "-"))
    tmpGenes <- as.character(unique(txDT$gene))
    tmpL <- split(txDT$gencoor, txDT$gene)[tmpGenes]
    tmpF <- parallel::mclapply(mc.cores = nCores, tmpGenes, function(gene_i){
        tmpO <- rep(NA, sum(txDT$gene == gene_i))
        txStr <- CDS_tab[CDS_tab$gene == gene_i,]$strand
        if(gene_i %in% CG$name){
            if(txStr == "+"){
                txSta <- CDS_tab[CDS_tab$gene == gene_i,]$start
                txEnd <- CDS_tab[CDS_tab$gene == gene_i,]$end
                tmpO[tmpL[[gene_i]] < txSta] <- "5'UTR"
                tmpO[tmpL[[gene_i]] >= txSta & tmpL[[gene_i]] <= txEnd] <- "CDS"
                tmpO[tmpL[[gene_i]] > txEnd] <- "3'UTR"
            }else if(txStr == "-"){
                txSta <- CDS_tab[CDS_tab$gene == gene_i,]$start
                txEnd <- CDS_tab[CDS_tab$gene == gene_i,]$end
                tmpO[tmpL[[gene_i]] > txSta] <- "5'UTR"
                tmpO[tmpL[[gene_i]] <= txSta & tmpL[[gene_i]] >= txEnd] <- "CDS"
                tmpO[tmpL[[gene_i]] < txEnd] <- "3'UTR"
            }
        }else{
            tmpO <- rep("non-coding", length(tmpO))
        }
        return(tmpO)
    }) %>% do.call(what = "c")
    txDT$geneRegion <- factor(tmpF, levels = c("5'UTR", "CDS", "3'UTR", "non-coding"))
    return(txDT)
}

#' Add exon number
#'
#' @param txDT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability is not available in Windows OS.
#'
#' @return data.table
#' @export
#'
tx_add_exonNumber <- function(txDT, geneAnnot, nCores = 1){
    # TODO: Check table integrity
    geneAnnot <- tx_cut_geneAnnotBytxDT(txDT, geneAnnot)
    EXONS <- exonGRanges(geneAnnot)
    tmpL <- split(txDT$gencoor, txDT$gene)[geneAnnot$name]
    txDT$exonNumber <- parallel::mclapply(mc.cores = nCores, as.character(unique(txDT$gene)), function(gene_i){
        strand_i <- as.character(GenomicRanges::strand(geneAnnot[geneAnnot$name == gene_i]))
        if(strand_i == "+"){
            tmpA <- which(tmpL[[gene_i]] %in% GenomicRanges::end(EXONS[[gene_i]])) -
                which(tmpL[[gene_i]] %in% GenomicRanges::start(EXONS[[gene_i]])) + 1
        }else if(strand_i == "-"){
            tmpA <- which(tmpL[[gene_i]] %in% GenomicRanges::start(EXONS[[gene_i]])) -
                which(tmpL[[gene_i]] %in% GenomicRanges::end(EXONS[[gene_i]])) + 1
        }
        rep(seq(EXONS[[gene_i]]), tmpA)
    }) %>% do.call(what = "c")
    txDT
}


#' Add exon place
#'
#' Adds a column specifying if the exon is the first, intermediate, or last exon
#'
#' @param txDT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability is not available in Windows OS.
#'
#' @return data.table
#' @export
tx_add_exonPlace <- function(txDT, geneAnnot, nCores = 1){
    if(!"exonNumber" %in% names(txDT)){
        cat("Exon number is not annotated in txDT \nAnnotating exon number\n")
        txDT <- tx_add_exonNumber(txDT, geneAnnot, nCores)
    } # TODO: Check table integrity. (Only if tx_add_exonNumber() is not used) (ELSE)
    tmpL <- split(txDT$exonNumber, txDT$gene)[as.character(unique(txDT$gene))]
    txDT$exonPlace <- parallel::mclapply(mc.cores = nCores, tmpL, function(x){
        if(all(x == 1)){
            return(factor(rep("single", length(x)),
                          levels = c("first", "intermediate", "last", "single")))
        }else{
            tmpF <- factor(rep("intermediate", length(x)), levels = c("first", "intermediate", "last", "single"))
            tmpF[x == 1] <- "first"
            tmpF[x == max(x)] <- "last"
            return(tmpF)
        }
    }) %>% do.call(what = "c")
    txDT
}

#' Add relative position in transcript
#'
#' Adds a numeric variable which represents the relative position of the
#' nucleotide along the transcript in a scale of 0 (start) to 1 (end).
#'
#' @param txDT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param round_dig integer. Decimal places for output to be round to.
#'
#' @return data.table
#' @export
tx_add_relTxPos <- function(txDT, round_dig = 3){
    # TODO: Check txDT integrity
    # TODO: Split only "txcoor"
    lapply(tx_split_DT(txDT), function(txDT){
        txDT$relTxPos <- round((txDT$txcoor)/max(txDT$txcoor), digits = round_dig)
        txDT
    }) %>% tx_merge_DT()
}

# Get functions ################################################################

#' Get data from a position and their neighboring positions
#'
#' Get the data from positions marked with TRUE in a logical variable as the
#' center and the values in its vicinity, delimited by a flank size.
#' Useful to aggregate data surrounding specific sites as motif locations.
#'
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param logi_col character. Name of column of logical class, which indicates queried sites
#' @param values_col character
#' @param upFlank numeric. Up-stream flank length
#' @param doFlank numeric. Down-stream flank length
#' @param addRowNames logical. Set to TRUE to add rownames in format "gene:txcoor".
#' Default is FALSE.
#'
#' @return matrix
#' @export
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
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param logi_col character. Name of column of logical class, which indicates queried sites
#' @param upFlank numeric. Up-stream flank length
#' @param doFlank numeric. Down-stream flank length
#' @param addNames logical. Set to TRUE to add names for each sequence in format
#' "gene:txcoor". Default is FALSE.
#'
#' @return character
#' @export
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
#' @param DT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#'
#' @return named numeric. Length of genes as per times they appear in DT
#' @export
tx_get_geneLengths <- function(DT){
    tmpTb <- table(as.character(DT$gene))
    tmpTb %>% as.numeric() %>% magrittr::set_names(names(tmpTb))
}

#' Get metagene at CDS
#'
#' @param txDT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().
#' @param colVars character. Names of columns for which values will be extracted
#' @param CDS_align character. Either "start", "end", or "spliceSite" depending on the desired
#' alignment, either to CDS start, CDS end, or splicing sites, respectively.
#' @param upFlank numeric. Up-stream flank length
#' @param doFlank numeric. Down-stream flank length
#'
#' @return list of matrices for each colVar
#' @export
tx_get_metageneAtCDS <- function(txDT, geneAnnot, colVars, CDS_align, upFlank, doFlank){
    check_GA_txDT_compat(txDT, geneAnnot)
    invisible(sapply(c("gencoor", "strand", "gene"), function(x) check_DThasCol(txDT, x)))
    geneAnnot <- geneAnnot[geneAnnot$name %in% unique(txDT$gene)]
    NCG <- geneAnnot[GenomicRanges::width(geneAnnot$thick) == 0]
    CG <- geneAnnot[GenomicRanges::width(geneAnnot$thick) > 0]
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
                                        addRowNames = TRUE)
        })
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
    }else if(CDS_align == "spliceSite"){
        txDT <- tx_add_spliceSitesLogical(txDT, geneAnnot)
        tmpFlanks <- lapply(colVars, function(colVar){
            tx_get_flanksFromLogicAnnot(DT = txDT,
                                        logi_col = "spliceSite",
                                        upFlank = upFlank,
                                        doFlank = doFlank,
                                        values_col = colVar,
                                        addRowNames = TRUE)})
    }else {stop("CDS_align should be either 'start' or 'end'.")}
    return(tmpFlanks %>% magrittr::set_names(colVars))
}

#' Get transcriptome sequences
#'
#' @param genome list. The full reference genome sequences, as loaded by
#' \code{\link{tx_load_genome}}() or prepackaged by BSgenome, see ?BSgenome::available.genomes
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().
#' @param outFile character. If specified, sequences will be written to this file name.
#' Otherwise a character vector will be output
#' @param nCores integer. Number of cores to run the function with. Multi-core
#' capability not available in Windows OS.
#'
#' @return character
#' @export
tx_get_transcriptSeqs <- function(genome, geneAnnot, outFile = NULL, nCores = 1){
    check_GA_genome_chrCompat(geneAnnot, genome)
    CHRS <- as.character(unique(GenomicRanges::seqnames(geneAnnot)))
    allSEQS <- parallel::mclapply(mc.cores = nCores, CHRS, function(iChr){
        subTXOME <- geneAnnot[which(as.logical(GenomicRanges::seqnames(geneAnnot) == iChr))]
        iBlocks <- IRanges::shift(S4Vectors::mcols(subTXOME)$blocks,
                                  IRanges::start(subTXOME) - 1)
        if(class(iBlocks)[1] != "CompressedIRangesList"){
            iBlocks <- IRanges::IRangesList(iBlocks, compress = TRUE) # Ensure CompressedIRangesList class
        }
        tmp <- stringr::str_sub(as.character(genome[[iChr]]), start = IRanges::start(unlist(iBlocks)),
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
        allSEQS[geneAnnot$name]
    }else{
        Biostrings::writeXStringSet(allSEQS[geneAnnot$name], filepath = outFile, format = "fasta")
    }
}

#' Get metagene regions
#'
#' Summarizes the values of a column(s) into a selected number of bins per
#' gene structure region (5'UTR, CDS, 3'UTR).
#'
#' @param txDT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().
#' @param colVars character. Names of columns for which values will be extracted
#' @param nBins_5UTR integer. Number of bins into which allocate data on 5'UTR regions
#' @param nBins_CDS integer. Number of bins into which allocate data on CDS regions
#' @param nBins_3UTR  integer. Number of bins into which allocate data on 3'UTR regions
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability not available in Windows OS.
#'
#' @return list of matrices for each colVar
#' @export
tx_get_metageneRegions <- function(txDT, geneAnnot, colVars, nBins_5UTR,
                                   nBins_CDS = NULL, nBins_3UTR = NULL,
                                   nCores = 1){
    warn_nCores(nCores)
    check_GA_txDT_compat(txDT, geneAnnot)
    invisible(sapply(c("gencoor", "strand", "gene"), function(x) check_DThasCol(txDT, x)))
    geneAnnot <- geneAnnot[geneAnnot$name %in% unique(txDT$gene)]
    NCG <- geneAnnot[GenomicRanges::width(geneAnnot$thick) == 0]
    CG <- geneAnnot[GenomicRanges::width(geneAnnot$thick) > 0]
    if(length(NCG)>0){warning(length(NCG), " non-coding genes where ommitted from analysis.")}
    if(sum(GenomicRanges::strand(CG) == "*") > 0){
        stop("Genes with no set strand (*) are not allowed in geneAnnot.")
    }
    if(is.null(nBins_CDS)){nBins_CDS <- nBins_5UTR * 9}
    if(is.null(nBins_3UTR)){nBins_3UTR <- nBins_5UTR * 9}
    genesDT <- unique(txDT$gene)
    if(!all(genesDT %in% geneAnnot$name)){
        stop("Gene annotation does not contain all genes in txDT.")
    }
    geneAnnot <- geneAnnot[geneAnnot$name %in% genesDT]
    if(!("geneRegion" %in% names(txDT))){
        cat("Gene regions are not annotated in txDT \nAnnotating gene regions\n")
        txDT <- tx_add_geneRegion(txDT, geneAnnot, nCores = nCores)
    }
    # Remove genes with small regions
    tmpR <- lapply(split(txDT$geneRegion, txDT$gene), FUN = "table") %>% do.call(what = "rbind") %>% data.frame()
    tmpR$gene <- rownames(tmpR)
    remGenes <- tmpR[tmpR$X3.UTR < nBins_3UTR | tmpR$X5.UTR < nBins_5UTR | tmpR$CDS < nBins_CDS, ]$gene
    txDT <- txDT[!(txDT$gene %in% remGenes), ]
    if(length(remGenes) > 0){
        warning(length(remGenes), " genes where removed from the analysis due to shorter length than selected binSizes.")
    }
    metageneBinNames <- c(paste("UTR5", seq(nBins_5UTR), sep = "_"),
                          paste("CDS", seq(nBins_CDS), sep = "_"),
                          paste("UTR3", seq(nBins_3UTR), sep = "_"))
    lapply(colVars, function(selCol){
        parallel::mclapply(mc.cores = nCores, unique(txDT$gene), function(gene_i){
            x <- split(txDT[txDT$gene %in% gene_i, ][[selCol]], txDT[txDT$gene %in% gene_i, "geneRegion"])
            tmpA <- c(tapply(x$`5'UTR`, cut(seq_along(x$`5'UTR`), breaks = nBins_5UTR),
                             function(a){mean(a, na.rm = TRUE)}),
                      tapply(x$CDS, cut(seq_along(x$CDS), breaks = nBins_CDS),
                             function(a){mean(a, na.rm = TRUE)}),
                      tapply(x$`3'UTR`, cut(seq_along(x$`3'UTR`), breaks = nBins_3UTR),
                             function(a){mean(a, na.rm = TRUE)}))
            tmpA[is.nan(tmpA)] <- NA
            tmpA
        }) %>% do.call(what = "rbind") %>% magrittr::set_colnames(metageneBinNames)
    }) %>% magrittr::set_names(colVars)
}

#' Get metagene by exons
#'
#' @param txDT data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
#' \code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.
#' @param colVars character. Names of columns for which values will be extracted
#' @param nBins integer. Number of bins into which exon data will be allocated.
#' @param geneAnnot GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().
#' @param rm_NArows logical. Remove rows of final matrix which consist of all NA values.
#' @param nCores integer. Number of cores to run the function with. Multicore
#' capability not available in Windows OS.
#'
#' @return list of matrices for each colVar
#' @export
tx_get_metageneExons <- function(txDT, colVars, nBins, geneAnnot = NULL, rm_NArows = TRUE, nCores = 1){
    if(!"exonNumber" %in% names(txDT)){
        warning("exonNumber is needed in txDT to extract exons, this function ",
                "calculates it but discards it everytime. If you want to keep ",
                "them please add it to the txDT by using tx_add_exonNumber()")
        if(is.null(geneAnnot)){stop("geneAnnot is needed to calculate exonNumber")}
        txDT <- tx_add_exonNumber(txDT = txDT, geneAnnot = geneAnnot, nCores = nCores)
    }
    metageneBinNames <- paste("exonBin", seq(nBins), sep = "_")
    lapply(colVars, function(selCol){
        tmpA <- parallel::mclapply(mc.cores = nCores, unique(txDT$gene), function(gene_i){
            x <- split(txDT[txDT$gene %in% gene_i, ][[selCol]], txDT[txDT$gene %in% gene_i, "exonNumber"])
            x <- x[sapply(x, "length") >= nBins]
            if(length(x) == 0){return(NULL)}
            lapply(x, function(x){
                tmpA <- tapply(x, cut(seq_along(x), breaks = nBins), FUN = function(a) mean(a, na.rm = TRUE))
                tmpA
            }) %>% do.call(what = "rbind") %>% magrittr::set_rownames(paste(gene_i, names(x), sep = "_"))
        }) %>% do.call(what = "rbind") %>% magrittr::set_colnames(metageneBinNames)
        tmpA[is.nan(tmpA)] <- NA
        if(rm_NArows){
            tmpA[rowMeans(is.na(tmpA)) != 1,]
        }else{
            tmpA
        }
    }) %>% magrittr::set_names(colVars)
}

# Statistical tests ############################################################

#' Likelihood Ratio Test
#'
#' @param DTL list. List of txDT, into which each element represents a replicate of an experimental dataset.
#' @param tVar character. Variable that represents number of trials, e.g. coverage
#' @param sVar character. Variable that represents number of successes, e.g. start_5p
#' @param test_groups factor. Factor specifying the "WT" and "KO" samples.
#' @param minTrials integer. Minimum number of trials in position to be considered in the statistical testing.
#'
#' @return data.table
#' @export
tx_test_LRTedgeR <- function(DTL, tVar, sVar, test_groups, minTrials = 50){
    DTL <- tx_unifyTxDTL(DTL)
    cMat_cov <- lapply(DTL, function(DT){
        DT[[tVar]]
    }) %>% do.call(what = "cbind") %>%
        magrittr::set_colnames(paste0(names(DTL), "-cov")) %>%
        magrittr::set_rownames(1:nrow(DTL[[1]]))
    cMat_sta <- lapply(DTL, function(DT){
        DT[[sVar]]
    }) %>% do.call(what = "cbind") %>%
        magrittr::set_colnames(paste0(names(DTL), "-sta")) %>% magrittr::set_rownames(1:nrow(DTL[[1]]))
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
    META <- data.table::data.table(names(DTL), group = test_groups)
    # design matrix
    designSL <- stats::model.matrix(~0+group, data = META)
    design <- edgeR::modelMatrixMeth(designSL)
    # Estimate dispersion
    cat("Estimating data dispersion and fitting GNB model\n")
    y1 <- edgeR::estimateDisp(y, design=design, trend = "none")
    # Fitting Negative Binomial Generalized Linear Models
    fit <- edgeR::glmFit(y1, design)
    #Contrast
    groupWT <- NULL #dummy variable
    groupKO <- NULL #dummy
    contr <- limma::makeContrasts(onlyIP = (groupWT-groupKO), levels = design)
    #Likelihood ratio test
    lrt <- edgeR::glmLRT(fit, contrast = contr[, "onlyIP"])
    #Results table
    RES_sta <- cbind(lrt$table, y$genes)
    RES_sta$FDR <- stats::p.adjust(RES_sta$PValue, method = "fdr")
    test_mat <- do.call(what = "cbind", lapply(DTL, function(DT) DT[[sVar]] / DT[[tVar]]))
    rownames(test_mat) <- paste(DTL[[1]][["gene"]], DTL[[1]][["txcoor"]], sep = ":")
    test_mat <- test_mat[paste(y$genes$gene, y$genes$txcoor, sep = ":"),] %>%
        rowMeansColG(test_groups)
    test_mat <- test_mat[,1] - test_mat[,2]
    RES_sta <- data.table::data.table(RES_sta, RD = test_mat)
    RES_sta <- RES_sta[order(RES_sta$FDR, decreasing = FALSE),]
    if("refSeq" %in% colnames(DTL[[1]])){
        RES_sta <- dplyr::right_join(RES_sta,
                                     dplyr::select(DTL[[1]],
                                                   c("chr", "gencoor", "strand",
                                                     "gene", "txcoor", "refSeq")),
                                     by = c("gene", "txcoor"))
    }else{
        RES_sta <- dplyr::right_join(RES_sta,
                                     dplyr::select(DTL[[1]],
                                                   c("chr", "gencoor", "strand",
                                                     "gene", "txcoor")))
    }
    return(RES_sta)
}


#' t-test in txDT list
#'
#' Apply a t-test per nucleotide position separating by groups over a txDT list.
#' At least 2 samples in each group are needed.
#'
#' @param DTL list. List of txDT, into which each element represents a replicate of an experimental dataset.
#' @param cont_var character. Name of column, specifying the continuous variable
#' to be used in test.
#' @param test_groups factor. Factor specifying the "WT" and "KO" samples.
#' @param test_na.rm logical. Remove NAs from tests
#' @param ... For use by \code{\link[genefilter]{rowttests}}
#'
#' @return data.table
#' @export
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
#' extends in leght of the second argument into both positive and negative
#' directions.
#'
#' @param position integer. Center of sequence
#' @param windowLength integer. Length of both downstream and upstream flanks
#'
#' @return integer
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
tx_dm3_geneAnnot <- function(){
    system.file("extdata", "toyGeneAnnot_Dmelan_chr4.bed", package = "txtools")
}

#' Load data - Case study #2
#'
#' @return character
#' @export
#'
tx_data_caseStudy2 <- function(){
    tmpF <- tempfile()
    utils::download.file(url = "https://drive.google.com/uc?export=download&id=1ICdtTNU0qK7ZetwNp5RzHDrYT5mrt6y2",
                         destfile = tmpF)
    readRDS(tmpF)
}

#' Download genome to temporary file - Case study # 1
#'
#' Sk1 genome. Source Schwartz et al., 2013 https://doi.org/10.1016/j.cell.2013.10.047
#'
#' @return character
#' @export
#'
sc_faGenome <- function(){
    tmpF <- tempfile()
    utils::download.file(url = "https://drive.google.com/uc?export=download&id=1IgUO5CKdHJh_2L2lp4ADfmGJlhCbKM9r",
                         destfile = tmpF)
    cat("sc_faGenome file downloaded to temporary file ", tmpF, "\n")
    tmpF
}

#' Download genome to temporary file - Case study # 3
#'
#' Thermococcus kodakarensis complete genome
#' Source: https://www.ncbi.nlm.nih.gov/assembly/GCF_000009965.1/
#'
#' @return character
#' @export
tk_faGenome <- function(){
    tmpF <- tempfile()
    utils::download.file(url = "https://drive.google.com/uc?export=download&id=1Iir-E0kVJ3RMHNsbOdSgi847t5lz1-re",
                         destfile = tmpF)
    cat("tk_faGenome file downloaded to temporary file ", tmpF, "\n")
    tmpF
}

#' Download gene annotation - Case study #1
#'
#' Download rRNA gene annotation of Saccharomyces cerevisae strain SK1 to temporary file
#'
#' @return character
#' @export
sc_geneAnnot <- function(){
    tmpF <- tempfile()
    utils::download.file(url = "https://drive.google.com/uc?export=download&id=1IlycxDiZoPcbpOsn7wRoXSx-v-79qrIW",
                         destfile = tmpF)
    cat("sc_geneAnnot file downloaded to temporary file ", tmpF, "\n")
    tmpF
}

#' Download gene annotation - Case study # 2
#'
#' Download mm9 gene annotation to temporary file.
#' Selected genes (toy example).
#'
#' Source: UCSC genes https://genome.ucsc.edu/cgi-bin/hgTables
#'
#' @return character
#' @export
mm_geneAnnot <- function(){
    tmpF <- tempfile()
    utils::download.file(url = "https://drive.google.com/uc?export=download&id=1IsZryJRYDtgevc8qXd_7pH-KqiNT8tSs",
                         destfile = tmpF)
    cat("mm_geneAnnot file downloaded to temporary file ", tmpF, "\n")
    tmpF
}

#' Download gene annotation - Case study # 3
#'
#' Download T. kodakarensis rRNA gene annotation to temporary file
#' Source:
#' @return character
#' @export
tk_geneAnnot <- function(){
    tmpF <- tempfile()
    utils::download.file(url = "https://drive.google.com/uc?export=download&id=1ImtH_ln_HXku2tAh2SaU-Io16zAV3ahQ",
                         destfile = tmpF)
    cat("tk_geneAnnot file downloaded to temporary file ", tmpF, "\n")
    tmpF
}
