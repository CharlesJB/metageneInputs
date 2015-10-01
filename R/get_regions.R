#' Extract the enhancers regions
#'
#' This uses the enhancers from the \code{FantomEnhancers.hg19} package. The
#' TPM values are calculated using the \code{FantomTSS.hg19} package by
#' calculating the mean expression of the TSS in each region.
#'
#' Finally, only the "robust" enhancers are returned, which means that only the
#' regions that overlaps at least one TSS from \code{FantomTSS.hg19}.
#'
#' @param cell_line The cell line to extract
#' @param keep_y Keep chrY. \code{TRUE} or \code{FALSE}. Default: \code{TRUE}.
#'
#' @return A \code{GRanges} corresponding to the enhancers
#'
#' @importFrom FantomEnhancers.hg19 get_fantom_enhancers_tpm
#' @import GenomeInfoDb
#' @import GenomicRanges
#'
#' @export
get_robust_enhancer_regions <- function(cell_line, keep_y = TRUE) {
  # Sanity check
  stopifnot(length(cell_line) == 1)
  stopifnot(is.logical(keep_y))

  # Prepare enhancers
  enhancers <-
    FantomEnhancers.hg19::get_fantom_enhancers_tpm(cell_lines = cell_line,
                                                   merge.FUN = mean)
  enhancers <- GenomicRanges::resize(enhancers, 2000, fix = "center")
  enhancers <- .clean_seqlevels(enhancers, keep_y)
  fetch_tpm(enhancers, cell_line)
}

.clean_seqlevels <- function(gr, keep_y) {
  gr <- GenomeInfoDb::keepStandardChromosomes(gr)
  gr <- GenomeInfoDb::sortSeqlevels(gr)
  if (keep_y == FALSE) {
    gr <- GenomeInfoDb::dropSeqlevels(gr, "chrY")
  }
  gr
}

#' Extract the enhancers regions
#'
#' This uses the enhanceeers from the \code{FantomEnhancers.hg19} package. The
#' TPM values are calculated using the \code{FantomTSS.hg19} package by
#' calculating the mean expression of the TSS in each region.
#'
#' Finally, only the "robust" enhancers are returned, which means that only the
#' regions that overlaps at least one TSS from \code{FantomTSS.hg19}.
#'
#' @param cell_line The cell line to extract
#' @param keep_y Keep chrY. \code{TRUE} or \code{FALSE}. Default: \code{TRUE}.
#'
#' @return A \code{GRanges} corresponding to the enhancers
#'
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import metagene
#'
#' @export
get_robust_promoter_regions <- function(cell_line, keep_y = TRUE) {
  # Sanity check
  stopifnot(length(cell_line) == 1)
  stopifnot(is.logical(keep_y))

  # Prepare enhancers
  data(promoters_hg19, package = "metagene")
  promoters <- .clean_seqlevels(promoters_hg19, keep_y)
  fetch_tpm(promoters, cell_line)
}

#' Fetch the robust TPM values.
#'
#' This function will first fetch the TPM values that overlap the regions
#' of interest from the robust TSS dataset (FantomTSS.hg19). When there is more
#' than one TSS for a genomic region, the mean value will be used. Only the
#' regions that overlaps at least one robust TSS will be returned.
#'
#' @param region A \code{GRanges} object corresponding to the regions of
#'               interest.
#' @param cell_line The cell line to extract the robust TSS TPM values.
#'
#' @return A \code{GRanges} corresponding to the subset of the \code{region}
#'         that overlaps at least on robust TSS with the mean robust TSS values
#'         for each range.
#'
#' @importFrom FantomTSS.hg19 get_fantom_tss_tpm
#' @import GenomicRanges
#' @import S4Vectors
#'
#' @export
fetch_tpm <- function(region, cell_line) {
  # Prepare TSS
  fantom_tss <- FantomTSS.hg19::get_fantom_tss_tpm(cell_lines = cell_line,
                                                   merge.FUN = mean)

  # Add TPM value
  expr <- GenomicRanges::mcols(fantom_tss)[[1]]
  GenomicRanges::mcols(region)[[cell_line]] <- 0
  hts <- GenomicRanges::findOverlaps(region, fantom_tss)
  qhts <- S4Vectors::queryHits(hts)
  shts <- S4Vectors::subjectHits(hts)
  expr <- sapply(split(expr[shts], qhts), mean)
  stopifnot(length(expr) == length(unique(qhts)))
  stopifnot(identical(names(expr), as.character(unique(qhts))))
  GenomicRanges::mcols(region)[[cell_line]][unique(qhts)] <- expr

  # Keep only regions that overlap Fantom's TSS
  GenomicRanges::subsetByOverlaps(region, fantom_tss)
}

#' Get all the robust TSS TPM values
#'
#' @param cell_line Limit the value returned to those of a single cell line. If
#'                  \code{NULL}, all cell lines are used. Default: \code{NULL}.
#' @param exclude Exclude a cell line. If \code{NULL}, all cell lines are used.
#'                Default: \code{NULL}
#'
#' @return A numeric vector of all the robut TSS from Fantom5
#'
#' @importFrom FantomTSS.hg19 get_fantom_library_name
#' @importFrom FantomTSS.hg19 get_fantom_tss_tpm
#' @import S4Vectors
#' @import IRanges
#'
#' @export
get_all_tss_tpm <- function(cell_line = NULL, exclude = NULL) {
  fantom_tpm <-
      GenomicRanges::mcols(FantomTSS.hg19::get_fantom_tss_tpm(cell_line))
  if (is.null(cell_line)) {
    fantom_tpm <- fantom_tpm[ grepl("CNhs", colnames(fantom_tpm))]
  } else {
    fantom_tpm <- lapply(fantom_tpm, S4Vectors::Rle)
  }
  if (!is.null(exclude)) {
    stopifnot(is.character(exclude))
    stopifnot(length(exclude) == 1)
    exclude <- FantomTSS::get_fantom_library_name(exclude)
    i <- ! (colnames(fantom_tpm) %in% exclude)
    fantom_tpm <- fantom_tpm[,i]
  }
  fantom_tpm <- lapply(fantom_tpm, round, digits = 2)
  fantom_tpm <- IRanges::RleList(fantom_tpm)
  fantom_tpm <- S4Vectors::unlist(fantom_tpm, use.names = FALSE)
  sort(fantom_tpm)
}

#' Split a GRanges into quantiles.
#'
#' @param gr The \code{GRanges} to extract quantile info from.
#' @param gr The tresholds for the quantiles.
#' @param cell_line The name of the column to use for splitting \code{gr}. If
#'                  \code{NULL}, the first metadata column is used. Default:
#'                  \code{NULL}.
#'
#' @return A \code{GRangesList} of the original \code{GRanges} splitted into
#'         quantiles.
#'
#' @import GenomicRanges
#'
#' @export
split_by_quantiles <- function(gr, quantiles, cell_line = NULL) {
  if (is.null(cell_line)) {
    cell_line <- colnames(GenomicRanges::mcols(gr))[1]
  }
  gr$quantile <- 0
  for (value in names(quantiles)) {
    j <- GenomicRanges::mcols(gr)[[cell_line]] <= quantiles[value]
    gr$quantile[j] <- value
  }
  j <- GenomicRanges::mcols(gr)[[cell_line]] == 0
  gr$quantile[j] <- "0%"
  gr$quantile <- as.numeric(as.character(gsub("%", "", gr$quantile)))
  split(gr, gr$quantile)
}
