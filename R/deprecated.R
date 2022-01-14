#' @export
tx_coverageDT <- function(x, geneAnnot, genome = NULL, nCores = 1){
    .Deprecated("tx_makeDT_coverage")
    tx_makeDT_coverage(x, geneAnnot, genome, nCores)
}

#' @export
tx_nucFreqDT <- function(x, geneAnnot, simplify_IUPAC = "splitForceInt", genome = NULL, nCores = 1){
    .Deprecated("tx_makeDT_nucFreq")
    tx_makeDT_nucFreq(x, geneAnnot, simplify_IUPAC, genome, nCores)
}

#' @export
tx_covNucFreqDT <- function(x, geneAnnot, simplify_IUPAC = "splitForceInt", genome = NULL, nCores = 1){
    .Deprecated("tx_makeDT_covNucFreq")
    tx_makeDT_covNucFreq(x, geneAnnot, simplify_IUPAC, genome, nCores)
}

#' @export
tx_filter_max_width <- function(x, thr, nCores = 1){
    .Deprecated("tx_filter_maxWidth")
    tx_filter_maxWidth(x, thr, nCores)
}

#' @export
tx_reads_mc <- function(reads, geneAnnot, minReads = 50, withSeq = FALSE, verbose = TRUE, nCores = 1){
    .Deprecated("tx_reads")
    tx_reads(reads, geneAnnot, minReads, withSeq, verbose, nCores)
}
