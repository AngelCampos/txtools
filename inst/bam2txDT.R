#!/usr/bin/env Rscript
# Title: Process bam files using txtools to create data.tables ########
# Author: Miguel Angel Garcia-Campos https://angelcampos.github.io/ ############
suppressPackageStartupMessages(require(optparse)) 
t0 <- Sys.time() # Start time

# Parsing Arguments ############################################################
option_list = list(
    make_option(c("-i", "--BAMfile"), type = "character", default = NULL,
                help="Input Paired-end reads alignment file in BAM format",
                metavar="character"),
    make_option(c("-o", "--outFile"), type = "character", default = "auto",
                help="Output file name. By default automatically sets it to the same path as BAMfile but with extension *.txDT.rds instead of *.bam",
                metavar="character"),
    make_option(c("-p", "--pairedEnd"), type = "logical", default = TRUE,
                help="Set to 'FALSE' for loading single-end reads alignments [default = %default]",
                metavar="logical"),
    make_option(c("-g", "--BEDfile_geneAnnot"), type = "character", default = NULL,
                help="Gene annotation file in BED-12 or BED-6 format",
                metavar="character"),
    make_option(c("-f", "--FASTAfile"), type = "character", default = NULL,
                help="Genome file in fasta format. If left empty reference sequence is not added [default = %default]",
                metavar="character"),
    make_option(c("-d", "--DT_datatype"), type = "character", default = NULL,
                help="Data.table type to output: 'cov' = coverage, 'covNuc' = Coverage & Nucleotide Frequency",
                metavar="character"),
    make_option(c("-r", "--removeLongerThan"), type = "integer", default = NA,
                help="Maximum length of reads to output [default = %default]",
                metavar="integer"),
    make_option(c("-m", "--minReadsGene"), type = "integer", default = 50,
                help="Minimum number of reads to report a gene in annotation [default = %default]",
                metavar="numeric"),
    make_option(c("-s", "--strandMode"), type = "integer", default = 1,
                help= paste0("Strand mode for loading BAM file. [default = %default]", 
                    "\n                1 = Direction is given by read 1 e.g. Directional Illumina", 
                    "\n                2 = Direction is given by read 2 (or just inverted in single-end files), e.g. Illumina TruSeq PE",
                    "\n                0 = Strand is unspecified, works as both directions."),
                metavar="numeric"),
    make_option(c("-n", "--nCores"), type = "integer", default = 2,
                help="Number of cores to be used [default= %default]",
                metavar="numeric"),
    make_option(c("-y", "--yieldSize"), type = "integer", default = 100000,
                help="Number of BAM records to be processed at a time [default = %default]",
                metavar="numeric"),
    make_option(c("-v", "--verbose"), type = "logical", default = TRUE,
                help="Set to 'FALSE' for showing task progress information [default = %default]",
                metavar="logical")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

BAMfile <- opt$BAMfile
outName <- opt$outFile
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

#Loading/Installing packages ###################################################
txtools_minVersion <- "0.0.7.4"
installLoad_CRAN <- function(package){
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
        install.packages(package, dependencies = TRUE, 
                         repos = "http://cran.us.r-project.org")
        library(package, character.only = TRUE, quietly = TRUE)
    }
}
suppressPackageStartupMessages(installLoad_CRAN("BiocManager"))
suppressPackageStartupMessages(installLoad_CRAN("magrittr"))
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
    if(!grepl(pattern = ".bam$", x = BAMfile)){stop("BAMfile, doesn't include extension .bam")}
    outName <- gsub(x = BAMfile, pattern = ".bam$", replacement = ".txDT.rds", perl = TRUE)
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
bam <- tx_load_bam(file = BAMfile,
                   yieldSize = ySize,
                   scanFlag = "default",
                   loadSeq = lSeq,
                   verbose = verb,
                   strandMode = strM,
                   pairedEnd = paired)
t1 <- Sys.time() # Loading BAM file time

# Convert to transcriptomic
txReads <- tx_reads(reads = bam,
                    geneAnnot = gA,
                    nCores = nCores,
                    minReads = minR,
                    withSeq = lSeq,
                    verbose = verb)

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
timeBam <- t1 - t0 # Total time taken
timePrc <- t2 - t1 # Total time taken
timeTot <- t2 - t0 # Total time taken
reportName <- strsplit(BAMfile, split = "/") %>% unlist %>% tail(1) %>% 
    gsub(x = ., pattern = ".bam$", replacement = ".txDT.log", perl = T)
readsInOut <- parallel::mclapply(mc.cores = nCores, txReads, names) %>% unlist
uniqReadsInOut <- unique(readsInOut)
report <- c("BAM file name:", BAMfile,
            "Paired-end reads in BAM file:", length(bam),
            "Output contains:", " ",
            "Number of genes:", length(unique(OUT$gene)),
            "Number of reads in output:", length(readsInOut),
            "Number of unique reads in output:", length(uniqReadsInOut),
            "Fraction of total reads in output:", round(length(uniqReadsInOut)/length(bam), 4),
            "Loading BAM time:", paste(round(timeBam, 2), units(timeBam), sep = " "),
            "Processing time:", paste(round(timePrc, 2), units(timePrc), sep = " "),
            "Total time taken:", paste(round(timeTot, 2), units(timeTot), sep = " ")) %>% 
    matrix(ncol = 2, byrow =T)
write.table(x = report, 
            file = reportName, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE,
            col.names = FALSE)

if(verb){
    print("Done!")
}
