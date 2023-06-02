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

#' Add number of nucleotide reads different to the reference genome
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' tx_add_diffNucToRef() was renamed to tx_add_misincCount()
#'
#' Add a column to DT of the sum of nucleotide frequency different to the
#' reference sequence counting deletions, without considering 'N's nor inserts
#' '.' into the calculation.
#'
#' @param DT data.table. A data.table object. See
#' \code{\link{tx_makeDT_coverage}}, \code{\link{tx_makeDT_nucFreq}} or
#' \code{\link{tx_makeDT_covNucFreq}} functions.
#'
#' @return data.table
#' @export
#' @keywords internal
tx_add_diffNucToRef <- function(DT){
    lifecycle::deprecate_warn("0.0.7", "tx_add_diffNucToRef()", "tx_add_misincCount()")
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("diffToRef")
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

#' Add different nucleotide reads to total ratio
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' tx_add_diffNucToRefRatio() was renamed to tx_add_misincRate()
#'
#' Add a column to DT of the ratio of different nucleotides to the total of
#' meaningfuk nucleotide reads not counting undetermined 'N' and inserts '.'.
#'
#' @param DT data.table
#' @param addDiffandTotalCols Set to TRUE to add counts of total nucleotides
#' read (nucTotal) and different to reference nucleotide (diffToRef) columns.
#'
#' @return data.table
#' @export
#' @keywords internal
#'
#' @seealso \code{\link{tx_add_diffNucToRef}} and \code{\link{tx_add_nucTotal}}
tx_add_diffNucToRefRatio <- function(DT, addDiffandTotalCols = FALSE){
    lifecycle::deprecate_warn("0.0.7", "tx_add_diffNucToRefRatio()", "tx_add_misincRate()")
    DT <- check_DT(DT) %>% hlp_removeColumnIfPresent("diffToRefRatio")
    tmpDT <- tx_add_diffNucToRef(DT) %>% tx_add_nucTotal()
    tmp <- round(tmpDT$diffToRef / tmpDT$nucTotal, 6)
    if(addDiffandTotalCols){DT <- tmpDT}
    tibble::add_column(DT, diffToRefRatio = tmp)
}

#' @export
tx_add_misincorpRateNucSpec <- function(DT, refNuc, misNuc, minNucReads = 20){
    lifecycle::deprecate_warn("0.0.7", "tx_add_misincorpRateNucSpec()", "tx_add_misincRateNucSpec()")
    tx_add_misincRateNucSpec(DT, refNuc, misNuc, minNucReads = minNucReads)
}
