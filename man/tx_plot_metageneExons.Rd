% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tx_plots.R
\name{tx_plot_metageneExons}
\alias{tx_plot_metageneExons}
\title{Plot metagene exons}
\usage{
tx_plot_metageneExons(
  txDT,
  colVars,
  nBins,
  geneAnnot = NULL,
  nCores = 1,
  summ_fun = "mean",
  smooth = TRUE,
  spar = 0.5,
  plot_type = "lineplot",
  xLabelJump = 10,
  normalize = FALSE
)
}
\arguments{
\item{txDT}{data.table. A table as output by the \code{\link{tx_makeDT_coverage}}(),
\code{\link{tx_makeDT_nucFreq}}() or \code{\link{tx_makeDT_covNucFreq}}() functions.}

\item{colVars}{character. Names of columns for which values will be extracted}

\item{nBins}{integer. Number of bins into which exon data will be allocated.}

\item{geneAnnot}{GRanges. Gene annotation as loaded by \code{\link{tx_load_bed}}().}

\item{nCores}{integer. Number of cores to run the function with.}

\item{summ_fun}{character. Summarizing function either "sum" or "mean". Default: "mean".}

\item{smooth}{logical. Set to FALSE for not smoothing line.}

\item{spar}{numeric. Smoothing parameter, typically (but not necessarily) in (0,1].}

\item{plot_type}{character. Type of plot to be output, either "lineplot" or "boxplot".}

\item{xLabelJump}{Number of bins to be skipped for label plotting}

\item{normalize}{logical. If set to TRUE, values are normalized so that the
area bellow the curve approximates to 1 for each variable.}
}
\value{
ggplot
}
\description{
Plot metagene exons
}
