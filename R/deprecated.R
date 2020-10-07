#' @export
tx_coverageDT <- function(x, geneAnnot, genome = NULL, nCores = 1){
    .Deprecated("tx_makeDT_coverage")
    tx_makeDT_coverage(x, geneAnnot, genome, nCores)
}

#' @export
tx_nucFreqDT <- function(x, geneAnnot, simplify_IUPAC = "splitForceInt", nCores = 1){
    .Deprecated("tx_makeDT_nucFreq")
    tx_makeDT_nucFreq(x, geneAnnot, simplify_IUPAC, nCores)
}

#' @export
tx_covNucFreqDT <- function(x, geneAnnot, simplify_IUPAC = "splitForceInt", nCores = 1){
    .Deprecated("tx_makeDT_covNucFreq")
    tx_makeDT_covNucFreq(x, geneAnnot, nCores)
}

#' @export
tx_filter_max_width <- function(x, thr, nCores = 1){
    .Deprecated("tx_filter_maxWidth")
    check_integer_arg(nCores, "nCores")
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
