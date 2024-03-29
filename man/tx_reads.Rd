% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tx_core.R
\name{tx_reads}
\alias{tx_reads}
\alias{tx_reads_mc}
\alias{tx_flushUnassigned}
\title{Transcriptomic reads convertion}
\usage{
tx_reads(
  reads,
  geneAnnot,
  minReads = 50,
  withSeq = FALSE,
  verbose = TRUE,
  ignore.strand = FALSE,
  nCores = 1
)
}
\arguments{
\item{reads}{GAlignments or GAlignmentPairs. Genomic alignments to be processed}

\item{geneAnnot}{GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}()}

\item{minReads}{integer. Minimum number of alignments required to overlap a gene}

\item{withSeq}{logical. Set to TRUE if sequence should be preserved; 'reads'
object should contain sequences.}

\item{verbose}{logical. Set to FALSE to show less information.}

\item{ignore.strand}{logical. Set to TRUE to allow alignments in a gene
ignoring the strand of the alignment. False by default.}

\item{nCores}{integer. Number of cores to run the function with. Multi-core
capability not available in Windows OS.}
}
\value{
GRanges
}
\description{
Assign aligned reads to their respective gene models and convert their
positions into their corresponding transcriptomic coordinate system.
It also stitches together paired-end aligned reads into a single sequence
in which dots '.' separate Read1 and Read2 by their corresponding insert.
}
\details{
To retrieve unassigned alignments use the function tx_get_unassignedAlignments()
}
