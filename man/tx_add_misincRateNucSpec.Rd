% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tx_core.R
\name{tx_add_misincRateNucSpec}
\alias{tx_add_misincRateNucSpec}
\alias{tx_add_misincorpRateNucSpec}
\title{Nucleotide specific misincorporation rate}
\usage{
tx_add_misincRateNucSpec(DT, refNuc, misNuc, minNucReads = 20)
}
\arguments{
\item{DT}{data.table. A table as output by the
\code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.}

\item{refNuc}{Reference nucleotide See \code{\link{tx_add_refSeqDT}}()}

\item{misNuc}{Misincorporated nucleotide}

\item{minNucReads}{Minimum number of nucleotides read needed to calculate
misincorporation rate}
}
\value{
data.table
}
\description{
Calculates the misincorporation rate from a \bold{ref}erence nucleotide to an
specific \bold{mis}incorporated nucleotide. For example, when looking for
positions in cytidines that were read as thymine in RNA after bisulphite treatment.
}
