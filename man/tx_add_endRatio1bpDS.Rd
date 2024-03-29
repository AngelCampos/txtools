% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tx_core.R
\name{tx_add_endRatio1bpDS}
\alias{tx_add_endRatio1bpDS}
\title{Add ends to coverage ratio 1bp down-stream}
\usage{
tx_add_endRatio1bpDS(DT, minCov = 50)
}
\arguments{
\item{DT}{data.table. A table as output by the \code{\link{tx_makeDT_coverage}}()
or \code{\link{tx_makeDT_covNucFreq}}() functions.}

\item{minCov}{numeric. Minimum coverage required to output ratio. If coverage
is less then an NA is output in that position.}
}
\value{
data.table
}
\description{
Add a column to DT of the read-ends to coverage ratio, shifted 1 base-pair
downstream. This means that the last measurement in any gene is always an NA
to account that there was no measurement to be shifted.
}
