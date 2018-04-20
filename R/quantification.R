#' Suitable for ATAC-seq or DNase-seq.
#'
#' Name sort the bams or this will take a long time.
#'
#' @param gr GRanges object with the regions to quantify.
#' @param bamlist Preferably a named list holding bam file paths, or a character vector.
#' @param nthreads Number of processors to use. Passed to Rsubread::featureCounts.
#' @return A matrix holding count data with axes regions x samples.
#' @export
quantifyCutsites <- function(gr, bamlist, nthreads = 8) {
  GenomicRanges::strand(gr) <- "*"

  saf <- data.frame(GeneID = as.character(gr), Chr = GenomicRanges::seqnames(gr),
                    Start = GenomicRanges::start(gr), End = GenomicRanges::end(gr),
                    Strand = GenomicRanges::strand(gr)) %>% set_rownames(NULL)

  if (is.null(names(bamlist))) {
    names(bamlist) <- lapply(bamlist, strsplit, "/") %>%
      lapply(unlist) %>% lapply(tail, 1) %>% unlist %>%
      make.unique(sep = "_")
  }

  cts <- Rsubread::featureCounts(bamlist, annot.ext = saf,
                                 nthreads = nthreads, isPairedEnd = FALSE, read2pos = 5,
                                 allowMultiOverlap = TRUE)

  rawcts <- cts$counts

  colnames(rawcts) <- names(bamlist)

  return(rawcts)
}

#' Suitable for NGS in general.
#'
#' Name sort the bams or this will take a long time.
#'
#' @param gr GRanges object with the regions to quantify.
#' @param bamlist Preferably a named list holding bam file paths, or a character vector.
#' @param nthreads Number of processors to use. Passed to Rsubread::featureCounts.
#' @return A matrix holding count data with axes regions x samples.
#' @export
quantifyReads <- function(gr, bamlist, nthreads = 8) {
  GenomicRanges::strand(gr) <- "*"

  saf <- data.frame(GeneID = as.character(gr), Chr = GenomicRanges::seqnames(gr),
                    Start = GenomicRanges::start(gr), End = GenomicRanges::end(gr),
                    Strand = GenomicRanges::strand(gr)) %>% set_rownames(NULL)

  if (is.null(names(bamlist))) {
    names(bamlist) <- lapply(bamlist, strsplit, "/") %>%
      lapply(unlist) %>% lapply(tail, 1) %>% unlist %>%
      make.unique(sep = "_")
  }

  #### NO 5' collapse here
  cts <- Rsubread::featureCounts(bamlist, annot.ext = saf,
                                 nthreads = nthreads, isPairedEnd = FALSE,
                                 allowMultiOverlap = TRUE)

  rawcts <- cts$counts

  colnames(rawcts) <- names(bamlist)

  return(rawcts)
}
