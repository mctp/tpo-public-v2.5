% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/triage.R
\name{.updateTableFilter}
\alias{.updateTableFilter}
\title{Updates table filter column based on triaging result.}
\usage{
.updateTableFilter(rep, tri, in_place = TRUE)
}
\arguments{
\item{rep}{variant table}

\item{tri}{triaging table}

\item{in_place}{update data.table in-place?}
}
\value{
rep
}
\description{
Note: for DB-backed tables update is always in-place and
expects non-subsetted dbVault
}
