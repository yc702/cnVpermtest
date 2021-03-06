% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tcgacnv_ampdel.R
\name{tcgacnv_ampdel}
\alias{tcgacnv_ampdel}
\title{A CNA permutation function for CNA amplification and deletion}
\usage{
tcgacnv_ampdel(cnvdatalist, segmentmean_name, segment_cutoff, permtime,
  core)
}
\arguments{
\item{cnvdatalist}{A list of CNA data with information about Chromosome, start and end position and segmentmean.}

\item{segmentmean_name}{The name of CNA segment mean column in the cnvdatalist}

\item{segment_cutoff}{Cut off value (positive or negative) of CNA segments in CNA data lists for amplification and deletion}

\item{permtime}{Permutation time}

\item{core}{Go parallel number of cores using}
}
\description{
This function allows you to do permutation of CNA segments for CNA amplification and deletion.
}
