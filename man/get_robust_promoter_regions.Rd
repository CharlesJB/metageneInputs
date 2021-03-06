% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/get_regions.R
\name{get_robust_promoter_regions}
\alias{get_robust_promoter_regions}
\title{Extract the enhancers regions}
\usage{
get_robust_promoter_regions(cell_line, keep_y = TRUE)
}
\arguments{
\item{cell_line}{The cell line to extract}

\item{keep_y}{Keep chrY. \code{TRUE} or \code{FALSE}. Default: \code{TRUE}.}
}
\value{
A \code{GRanges} corresponding to the enhancers
}
\description{
This uses the enhanceeers from the \code{FantomEnhancers.hg19} package. The
TPM values are calculated using the \code{FantomTSS.hg19} package by
calculating the mean expression of the TSS in each region.
}
\details{
Finally, only the "robust" enhancers are returned, which means that only the
regions that overlaps at least one TSS from \code{FantomTSS.hg19}.
}

