% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mut_rain.R
\name{mc.mutRain}
\alias{mc.mutRain}
\title{Compute distance between mutations for rainfall plots}
\usage{
mc.mutRain(mutnet.obj, col.chrom = "chrom", col.pos = "pos",
  col.ref = "ref.allele", col.alt = "alt.allele", symetrize = c(TRUE,
  FALSE)[1])
}
\arguments{
\item{mutnet.obj}{A dataframe of mutation info}

\item{col.chrom}{Name of the column with the chromosomes}

\item{col.pos}{Name of the column with the mutations positions}

\item{col.ref}{Name of the column with the reference allele}

\item{col.alt}{Name of the column with the alternative allele}

\item{symetrize}{Boolean to symetrize reversable mutations or not}
}
\value{
A table with relative distance for each mutation
}
\description{
Compute distance between mutations for rainfall plots
}

