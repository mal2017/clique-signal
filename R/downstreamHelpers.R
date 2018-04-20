#' index_by_clique
#' @export
index_by_clique <- function(object) {
  stopifnot(class(object) == "CRCExperiment")
  cliquewise_tfbs <- tfbs_by_clique(object)
  # https://support.bioconductor.org/p/82452/
  message("Indexing accessible sites by clique-grouped TFBS...", appendLF = F)
  tictoc::tic()
  hits <- findOverlaps(rowRanges(object), cliquewise_tfbs)
  mat <- matrix(FALSE, nrow = S4Vectors::queryLength(hits),
                ncol = S4Vectors::subjectLength(hits))
  mat[S4Vectors::as.matrix(hits)] <- TRUE
  colnames(mat) <- names(cliquewise_tfbs)
  rownames(mat) <- rownames(object)
  tictoc::toc()
  Matrix::Matrix(mat, sparse=T)
}

#' get_diff_cliques
#' @export
get_diff_cliques <- function(object, ix, contrast) {
  # spec contrast of condition levels
  NULL
  ## set attributes so plotting function will take it
}


#' plot_cliques
#' @export
plot_cliques <- function(res, plot=c("volcano","waterfall","tsne","pca","clustergram"), interactive=T) {
  NULL
}
