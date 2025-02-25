% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vault_annotate.R
\name{getGenes}
\alias{getGenes}
\title{get all available genes for a sample}
\usage{
getGenes(SOM, GRM, STR, FUS, GENE.COPY)
}
\arguments{
\item{SOM}{annotated somatic}

\item{GRM}{annotated germline}

\item{STR}{annotated structural}

\item{FUS}{annotated fusions}

\item{GENE.COPY}{annotated gene copy number}
}
\value{
DF ensembl/symbol
}
\description{
get all available genes for a sample
}
