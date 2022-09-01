#' Retrieve dumped alignments
#'
#' Alias for GenomicAlignments::getDumpedAlignments. Retrieve ambiguous reads
#' not included in the loaded object. For more info check:
#' \code{\link[GenomicAlignments]{findMateAlignment}}
#'
#' @export
getDumpedAlignments <- function(){
    GenomicAlignments::getDumpedAlignments()
}

# Pipe operator
`%>%` <- magrittr::`%>%`
