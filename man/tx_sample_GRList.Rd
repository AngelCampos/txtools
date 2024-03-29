% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tx_core.R
\name{tx_sample_GRList}
\alias{tx_sample_GRList}
\title{Sampling alignments}
\usage{
tx_sample_GRList(x, p, nCores = 1)
}
\arguments{
\item{x}{CompressedGRangesList. A list containing alignments, meant to be used
for the output of the \code{\link{tx_reads}} function.}

\item{p}{Probabilty for each read to be sampled.}

\item{nCores}{integer. Number of cores to run the function with. Multicore
capability not available in Windows OS.}
}
\value{
CompressedGRangesList.
}
\description{
Sampling of alignmnets in a GRanges list using a binomial distribution.
}
