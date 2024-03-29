% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tx_plots.R
\name{tx_plot_metageneRegions}
\alias{tx_plot_metageneRegions}
\title{Plot metagene by regions}
\usage{
tx_plot_metageneRegions(
  txDT,
  geneAnnot,
  colVars,
  nBins_5UTR,
  nBins_CDS = NULL,
  nBins_3UTR = NULL,
  summ_fun = "mean",
  smooth = TRUE,
  spar = 0.5,
  nCores = 1,
  plot_type = "lineplot",
  normalize = FALSE
)
}
\arguments{
\item{txDT}{data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
\code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.}

\item{geneAnnot}{GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().}

\item{colVars}{character. Names of columns for which values will be extracted}

\item{nBins_5UTR}{integer. Number of bins into which allocate data on 5'UTR regions}

\item{nBins_CDS}{integer. Number of bins into which allocate data on CDS regions}

\item{nBins_3UTR}{integer. Number of bins into which allocate data on 3'UTR regions}

\item{summ_fun}{character. Summarizing function either "sum" or "mean". Default: "mean"}

\item{smooth}{logical. Set to FALSE for not smoothing line.}

\item{spar}{numeric. Smoothing parameter, typically (but not necessarily) in (0,1].}

\item{nCores}{integer. Number of cores to run the function with. Multicore
capability not available in Windows OS. @param plot_type character. Type of plot to be output, either "lineplot" or
"boxplot".}

\item{plot_type}{character. Type of plot to be output, either "lineplot" or "boxplot".}

\item{normalize}{logical. If set to TRUE, values are normalized so that the
area bellow the curve approximates to 1 for each variable.}
}
\value{
ggplot
}
\description{
Plot metagene by regions
}
