#' Get the available targets for a cell_line
#'
#' This function will extract the names of all the targets that have at least
#' one bam file available on the ENCODE web site.
#'
#' @param cell_line The name of the cell line (biosample).
#'
#' @return A vector containing the names of the targets
#'
#' @import ENCODExplorer
#'
#' @export
get_encode_targets <- function(cell_line) {
    qry <- ENCODExplorer::queryEncode(assay = "ChIP-seq", biosample = cell_line,
                                      organism = "Homo sapiens",
                                      file_format = "bam")
    targets <- unique(vapply(strsplit(qry[["experiment"]][["target"]], "-"),
           `[`, 1, FUN.VALUE = character(1)))
    targets[order(targets)]
}

#' Produce the design for a specific target
#'
#' @param target The name of the target.
#'
#' @return The design in \code{data.frame} format.
#'
#' @import ENCODExplorer
#'
#' @export
get_encode_design <- function(target, cell_line) {
    chip <- ENCODExplorer::queryEncode(assay = "ChIP-seq",
                                       biosample = cell_line,
                                       file_format = "bam",
                                       target = paste0(target, "-human"),
                                       status = "released")
    ctrl_accession <- vapply(strsplit(chip[["experiment"]][["controls"]], "/"),
                             `[`, 3, FUN.VALUE = character(1))
    ctrl_accession <- unique(ctrl_accession)
    ctrl <- ENCODExplorer::queryEncode(assay = "ChIP-seq",
                                       biosample = cell_line,
                                       file_format = "bam",
                                       status = "released")
    i <- ctrl[["experiment"]][["accession"]] %in% ctrl_accession
    ctrl[["experiment"]] <- ctrl[["experiment"]][i,]
    chip_bam <- paste0(chip[["experiment"]][["file_accession"]], ".bam")
    ctrl_bam <- paste0(ctrl[["experiment"]][["file_accession"]], ".bam")
    bam_files <- c(chip_bam, ctrl_bam)
    i <- bam_files %in% chip_bam
    j <- bam_files %in% ctrl_bam
    design <- data.frame(bam_files = bam_files)
    design[[target]] <- 0
    design[[target]][i] <- 1
    design[[target]][j] <- 2
    stopifnot(all(design[[target]] != 0))
    design
}
