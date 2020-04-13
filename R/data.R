#' IUPAC code extended to consider deletions
#'
#' The IUPAC_CODE_MAP_extended named character vector contains the mapping from
#' the IUPAC nucleotide ambiguity codes to their meaning adding the deletion
#' code "-", this is required for computing the consensus strings between
#' overlapping reads.
#'
#' @format Named character vector
#'
#' @source Adapted from \code{\link[Biostrings]{IUPAC_CODE_MAP}}
"IUPAC_CODE_MAP_extended"

#' IUPAC ambiguity alphabet (2nc)
#'
#' IUPAC nucleotide ambiguity alphabet for 2 combinations of two nucleotides,
#' for use mostly by the txtools functions that generate datatables with
#' nucleotide frequency
#'
#' @format Character vector
#'
"IUPAC_code_2nucs"

#' Simplified nucleotide alphabet
#'
#' Simplified nucleotide alphabet for use mostly by the txtools functions that
#' generate datatables with nucleotide frequency
#'
#' @format Character vector
#'
"IUPAC_code_simpl"
