#' Quantify reads in ATAC-seq or DNase-seq bams.
#'
#' @param gr GRanges object with the regions to quantify.
#' @param bamlist Preferably a named list holding bam file paths, or a character vector.
#' @param nthreads Number of processors to use. Passed to Rsubread::featureCounts.
#' @return A matrix holding count data with axes regions x samples.
#' @export
quantifyCutsites <- function(gr, bamlist, nthreads = 8) {
    GenomicRanges::strand(gr) <- "*"

    saf <- data.frame(GeneID = as.character(gr), Chr = GenomicRanges::seqnames(gr), Start = GenomicRanges::start(gr),
        End = GenomicRanges::end(gr), Strand = GenomicRanges::strand(gr)) %>% set_rownames(NULL)

    if (is.null(names(bamlist))) {
        names(bamlist) <- lapply(bamlist, strsplit, "/") %>% lapply(unlist) %>% lapply(tail, 1) %>%
            unlist %>% make.unique(sep = "_")
    }

    cts <- Rsubread::featureCounts(bamlist, annot.ext = saf, nthreads = nthreads, isPairedEnd = FALSE,
        read2pos = 5, allowMultiOverlap = TRUE)

    aligned_depth <- cts$stat %>%
      subset(Status != "Unassigned_Unmapped") %>%
      .[,-c(1)] %>% lapply(sum)

    rawcts <- cts$counts

    colnames(rawcts) <- names(bamlist)
    names(aligned_depth) <- names(bamlist)

    return(list(cts = rawcts, aligned_depth = aligned_depth))
}

#' Suitable for NGS in general.
#'
#' @param gr GRanges object with the regions to quantify.
#' @param bamlist Preferably a named list holding bam file paths, or a character vector.
#' @param nthreads Number of processors to use. Passed to Rsubread::featureCounts.
#' @param paired_end Treat bams as paired-end or not. Defaults to TRUE.
#' @return A matrix holding count data with axes regions x samples.
#' @export
quantifyReads <- function(gr, bamlist, nthreads = 8, paired_end = T) {
    GenomicRanges::strand(gr) <- "*"

    saf <- data.frame(GeneID = as.character(gr), Chr = GenomicRanges::seqnames(gr), Start = GenomicRanges::start(gr),
        End = GenomicRanges::end(gr), Strand = GenomicRanges::strand(gr)) %>% set_rownames(NULL)

    if (is.null(names(bamlist))) {
        names(bamlist) <- lapply(bamlist, strsplit, "/") %>% lapply(unlist) %>% lapply(tail, 1) %>%
            unlist %>% make.unique(sep = "_")
    }

    #### NO 5' collapse here
    cts <- Rsubread::featureCounts(bamlist, annot.ext = saf, nthreads = nthreads, isPairedEnd = paired_end,
        allowMultiOverlap = F, largestOverlap = T, requireBothEndsMapped = F)

    aligned_depth <- cts$stat %>%
      subset(Status != "Unassigned_Unmapped") %>%
      .[,-c(1)] %>% lapply(sum)

    rawcts <- cts$counts

    colnames(rawcts) <- names(bamlist)
    names(aligned_depth) <- names(bamlist)

    return(list(cts = rawcts, aligned_depth = aligned_depth))
}
