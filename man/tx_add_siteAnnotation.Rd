% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tx_core.R
\name{tx_add_siteAnnotation}
\alias{tx_add_siteAnnotation}
\title{Add 1bp-site logical annotation}
\usage{
tx_add_siteAnnotation(DT, GRanges, colName, nCores = 1)
}
\arguments{
\item{DT}{data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
\code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.}

\item{GRanges}{GRanges. Ranges of length 1, to be marked in the data.table}

\item{colName}{character. Name of the new column to be added.}

\item{nCores}{integer. Number of cores to run the function with. Multicore
capability not available in Windows OS.}
}
\value{
data.table
}
\description{
Add a logical variable column in DT, for which the genomic coordinates of a
GRanges object are used to annotate 1 bp sites.
In which TRUE equals to presence of the site in the provided GRanges object.
This annotation is useful, for example, when marking already known RNA
modification sites.
}
