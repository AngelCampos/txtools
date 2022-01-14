#' IUPAC code extended to consider deletions
#'
#' The IUPAC_CODE_MAP_extended named character vector contains the mapping from
#' the IUPAC nucleotide ambiguity codes to their meaning adding the deletion
#' code "-", this is required for computing the consensus strings between
#' overlapping reads.
#'
#' @format Named character vector
#'
#' @source Adapted from Biostrings object \code{\link[Biostrings]{IUPAC_CODE_MAP}}
"IUPAC_CODE_MAP_extended"

#' IUPAC ambiguity alphabet (2nc)
#'
#' IUPAC nucleotide ambiguity alphabet for combinations of two nucleotides,
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

#' Genomic Alignments example - Sk1 yeast
#'
#' Paired-end alignment data extracted from the Schwartz et al., 2013, Sk1 yeast strain dataset.
#'
#' @format GAlignmentPairs
#'
"bam_sk1"


#' Gene annotation example - Sk1 yeast
#'
#' 2 genes models for the yeast strain Sk1
#'
#' @format GRanges
#'
"gA_sk1"

#' Genome example - Sk1 yeast
#'
#' Genomic sequence for chromosomes 4 and 5 of the yeast strain Sk1
#'
#' @format DNAStringSet
#'
"genome_sk1"
