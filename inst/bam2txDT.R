#!/usr/bin/env Rscript
# Title: Process bam files using txtools to create data.tables ########
# Author: Miguel Angel Garcia-Campos https://angelcampos.github.io/ ############
suppressPackageStartupMessages(require(optparse)) 
t0 <- Sys.time() # Start time

# Parsing Arguments ############################################################
option_list = list(
    make_option(c("-i", "--BAMfile"), type = "character", default = NULL,
                help="Input mapped/algined reads file in BAM format",
                metavar="character"),
    make_option(c("-p", "--pairedEnd"), type = "logical", default = TRUE,
                help="Set to 'FALSE' for loading single-end reads alignments [default = %default]",
                metavar="logical"),
    make_option(c("-g", "--BEDfile_geneAnnot"), type = "character", default = NULL,
                help="Gene annotation file in BED-12 or BED-6 format (no introns)",
                metavar="character"),
    make_option(c("-f", "--FASTAfile"), type = "character", default = NULL,
                help="Genome file in fasta format. If left empty the reference sequence will not be added [default = %default]",
                metavar="character"),
    make_option(c("-d", "--DT_datatype"), type = "character", default = NULL,
                help="Data.table type to output: 'cov' = coverage, 'covNuc' = Coverage & Nucleotide Frequency",
                metavar="character"),
    make_option(c("-o", "--outName"), type = "character", default = "auto",
                help="Output name, expected to end in '.rds'. By default the BAM file name will be used changing '.bam' for '.txDT.rds'",
                metavar="character"),
    make_option(c("-r", "--removeLongerThan"), type = "integer", default = NA,
                help="Maximum length of reads to output [default = %default]. By default will not remove reads.",
                metavar="integer"),
    make_option(c("-m", "--minReadsGene"), type = "integer", default = 50,
                help="Minimum number of reads to report a gene in annotation [default = %default]",
                metavar="integer"),
    make_option(c("-R", "--rescueSingletons"), type = "logical", default = FALSE,
                help="Separate mapped singletons and process separately to rescue those mapped reads (only for pairedEnd == TRUE). [default = %default]",
                metavar="logical"),
    make_option(c("-s", "--strandMode"), type = "integer", default = 1,
                help= paste0("Strand mode for loading BAM file. [default = %default]", 
                             "\n                1 = Direction is given by read 1 e.g. Directional Illumina", 
                             "\n                2 = Direction is given by read 2 (or just inverted in single-end files), e.g. Illumina TruSeq PE"),
                metavar="integer"),
    make_option(c("-S", "--ignoreStrand"), type = "logical", default = FALSE,
                help="Ignore strand of mapping, use for unstranded library preparations. [default = %default]",
                metavar="logical"),
    make_option(c("-n", "--nCores"), type = "integer", default = 2,
                help="Number of cores to be used [default= %default]",
                metavar="integer"),
    make_option(c("-y", "--yieldSize"), type = "integer", default = 100000,
                help="Number of BAM records to be processed at a time [default = %default]",
                metavar="integer"),
    make_option(c("-v", "--verbose"), type = "logical", default = TRUE,
                help="Set to 'FALSE' for showing task progress information [default = %default]",
                metavar="logical")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

BAMfile <- opt$BAMfile
paired <- opt$pairedEnd
geneAnnot <- opt$BEDfile_geneAnnot
genome <- opt$FASTAfile
dtType <- opt$DT_datatype
remL <- opt$removeLongerThan
minR <- opt$minReadsGene
strM <- opt$strandMode
nCores <- opt$nCores
ySize <- opt$yieldSize
verb <- opt$verbose
rescueSingle <- opt$rescueSingletons
ignStrand <- opt$ignoreStrand
outName <- opt$outName

#Loading/Installing packages ###################################################
txtools_minVersion <- "0.1.2"
installLoad_CRAN <- function(package){
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
        install.packages(package, dependencies = TRUE, 
                         repos = "http://cran.us.r-project.org")
        library(package, character.only = TRUE, quietly = TRUE)
    }
}
installLoad_bioC <- function(package){
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
        BiocManager::install(package)
        library(package, character.only = TRUE, quietly = TRUE)
    }
}
suppressPackageStartupMessages(installLoad_CRAN("BiocManager"))
suppressPackageStartupMessages(installLoad_CRAN("magrittr"))
suppressPackageStartupMessages(installLoad_CRAN("parallel"))
suppressPackageStartupMessages(installLoad_bioC("Rsamtools"))
suppressPackageStartupMessages(installLoad_bioC("GenomicRanges"))

if(!requireNamespace("txtools", quietly = TRUE)){
    BiocManager::install("AngelCampos/txtools")
}else if(packageVersion("txtools") < txtools_minVersion){
    BiocManager::install("AngelCampos/txtools")
    suppressPackageStartupMessages(library("txtools"))
}else{
    suppressPackageStartupMessages(library("txtools"))
}

# Process arguments ############################################################
# Load sequence or not
if(dtType == "cov"){
    lSeq <- FALSE
}else if(dtType == "covNuc"){
    lSeq <- TRUE
}else{
    stop("-d --DT_datatype argument must be one of the options: 'cov' or 'covNuc'")
}
# Check yieldSize to be integer
txtools:::check_integer_arg(ySize, "ySize")
# Set OUTPUT filename
if(outName == "auto"){
    outName <- basename(BAMfile) %>% 
        gsub(x = ., pattern = ".bam$", replacement = ".txDT.rds", perl = T)
}

# Output all genes even with no reads overlapping
if(minR == 0){makeFULL <- TRUE}else{makeFULL <- FALSE}
if(minR == 0){minR <- 1}

# MAIN program #################################################################
# Load gene annotation
gA <- tx_load_bed(geneAnnot)
if(verb){
    cat("Gene annotation loaded with", length(gA), "gene models.")
}
# Load genome or NULL
if(!is.null(genome)){
    GENOME <- tx_load_genome(genome)
}else{
    GENOME <- NULL
}

# Load BAM file
if(!rescueSingle){
    bam <- tx_load_bam(file = BAMfile,
                       yieldSize = ySize,
                       scanFlag = "default",
                       loadSeq = lSeq,
                       verbose = verb,
                       strandMode = strM,
                       pairedEnd = paired)
    t1 <- Sys.time() # Loading BAM file time
    bamLen <- length(bam)
}else if(rescueSingle){
    if(!paired){stop("--pairedEnd argument is expected to be TRUE if --rescueSingletons is set to TRUE\n")}
    filtParam_1 <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE, hasUnmappedMate = FALSE)) # Fully paired
    filtParam_2 <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE, hasUnmappedMate = TRUE)) # Rescuable
    tmpFile1 <- tempfile(); tmpFile2 <- tempfile()
    invisible(filterBam(file = BAMfile, destination = tmpFile1, param = filtParam_1))
    invisible(filterBam(file = BAMfile, destination = tmpFile2, param = filtParam_2))
    tmpBAM_1 <- tx_load_bam(file = tmpFile1, 
                            pairedEnd = TRUE,
                            loadSeq = lSeq,
                            verbose = verb, 
                            yieldSize = ySize, 
                            scanFlag = "default", 
                            strandMode = strM)
    if(verb){cat("\nLoading singletons bam file\n")}
    tmpBAM_2 <- tx_load_bam(file = tmpFile2,
                            pairedEnd = FALSE,
                            loadSeq = lSeq,
                            verbose = verb,
                            yieldSize = ySize, 
                            scanFlag = "default",
                            strandMode = strM)
    t1 <- Sys.time() # Loading BAM files
    invisible(file.remove(c(tmpFile1, tmpFile2)))
    bamLen <- length(tmpBAM_1) + length(tmpBAM_2)
}
# txtools:::hlpr_ReadsInGene
# Convert to transcriptomic
if(!rescueSingle){
    txReads <- tx_reads(reads = bam,
                        geneAnnot = gA,
                        nCores = nCores,
                        minReads = minR,
                        withSeq = lSeq,
                        verbose = verb, 
                        ignore.strand = ignStrand)
}else if(rescueSingle){
    tmp_txReads_1 <- tx_reads(reads = tmpBAM_1,
                              minReads = minR,
                              geneAnnot = gA,
                              nCores = nCores, 
                              withSeq = lSeq,
                              verbose = verb, 
                              ignore.strand = ignStrand)
    tmp_txReads_2 <- tx_reads(reads = tmpBAM_2, 
                              geneAnnot = gA, 
                              nCores = nCores, 
                              withSeq = lSeq,
                              verbose = verb, 
                              ignore.strand = ignStrand)
    # Merge lists gene-wise
    allGenes <- union(names(tmp_txReads_1), names(tmp_txReads_2))
    txReads <- parallel::mclapply(mc.cores = nCores, allGenes, function(gene){
        tmp1 <- tmp_txReads_1[[gene]]
        tmp2 <- tmp_txReads_2[[gene]]
        if(is.null(tmp1)){
            seqlevels(tmp2) <- allGenes
            return(tmp2)
        }else if(is.null(tmp2)){
            seqlevels(tmp1) <- allGenes
            return(tmp1)
        }else{
            seqlevels(tmp1) <- allGenes; seqlevels(tmp2) <- allGenes
            return(c(tmp1, tmp2))
        }
    }) %>% GenomicRanges::GRangesList() %>% 
        magrittr::set_names(allGenes) %>% 
        GenomicRanges::GRangesList()
}

# Filter by length
if(!is.na(remL)){
    txReads <- tx_filter_maxWidth(x = txReads, thr = remL, nCores = nCores)
}

# Data.table generation
if(dtType == "cov"){
    if(verb){
        cat("Generating coverage data.table")
    }
    OUT <- tx_makeDT_coverage(x = txReads, 
                              geneAnnot = gA, 
                              nCores = nCores,
                              fullDT = makeFULL,
                              genome = GENOME)
}else if(dtType == "covNuc"){
    if(verb){
        cat("Generating coverage and nucleotide frequency data.table. \n")
    }
    OUT <- tx_makeDT_covNucFreq(x = txReads, 
                                geneAnnot = gA,
                                nCores = nCores,
                                fullDT = makeFULL,
                                genome = GENOME)
}else{print("This message should not be printed, let the maintainer know.")}
t2 <- Sys.time() # Creating data.table time
# NOTE: In practice using too many cores made for longer processing times
# newNCores <- min(nCores, 10)

# Saving file as .rds
saveRDS(object = OUT, file = outName)

# Report #######################################################################
timeBam <- t1 - t0 # Time loading BAM
timePrc <- t2 - t1 # Time processing
timeTot <- t2 - t0 # Total time 

reportName <- basename(outName) %>% gsub(x = ., pattern = ".rds$", replacement = ".log", perl = T)

readsInOut <- parallel::mclapply(mc.cores = nCores, txReads, names) %>% unlist()
uniqReadsInOut <- unique(readsInOut)
report <- c("BAM file name:", BAMfile,
            "Reads in BAM file:", bamLen,
            "Output contains:", " ",
            "    Number of genes:", length(unique(OUT$gene)),
            "    Number of reads in output:", length(readsInOut),
            "    Number of unique reads in output:", length(uniqReadsInOut),
            "    Fraction of total reads in output:", round(length(uniqReadsInOut)/bamLen, 4),
            "Total time taken:", paste(round(timeTot, 2), units(timeTot), sep = " "),
            "    Loading BAM time:", paste(round(timeBam, 2), units(timeBam), sep = " "),
            "    Processing time:", paste(round(timePrc, 2), units(timePrc), sep = " ")) %>% 
    matrix(ncol = 2, byrow =T)
write.table(x = report, 
            file = reportName, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE,
            col.names = FALSE)

if(verb){print("bam2txDT.R finished!")}
