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
#' Supply a CRCExperiment, and index matrix. Other params passed
#' to chromVAR::differentialDeviations()
#
#' @export
get_diff_cliques <- function(object, ix, contrast = NULL, parametric = TRUE,
                             score_type = c("dev","z"),
                             return = c("tibble","all")) {
  score_type <- score_type[1]
  return <- return[1]
  stopifnot(nrow(object) == nrow(ix))
  stopifnot(length(contrast) == 2 | is.null(contrast))

  if (!is.null(contrast)) {
    object <- object[, object$CONDITION %in% contrast]
  }

  colData(object) <- droplevels(colData(object))
  message("Computing deviations via chromVAR... ", appendLF = F)
  tictoc::tic()
  dev <- chromVAR::computeDeviations(object, annotations = ix)
  diff_dev <- chromVAR::differentialDeviations(dev,groups = "CONDITION",
                                               alternative = "two.sided",
                                               parametric = parametric)
  tictoc::toc()

  if (score_type == "dev") {
    scores <- chromVAR::deviations(dev)
    scores <- tibble::as_tibble(scores,rownames="clique")
  } else if (score_type == "z") {
    scores <- chromVAR::deviationScores(dev)
    scores <- tibble::as_tibble(scores,rownames="clique")
  }

  split(rownames(colData(object)),
        colData(object)$CONDITION) -> names_by_cond

  diffs <- tibble::tibble(clique=scores$clique)
  for (i in 1:length(names_by_cond)) {
    grp_name <- names(names_by_cond)[i]
    samples <- names_by_cond[[i]]
    diffs[[grp_name]] <- rowMeans(scores[samples])
  }

  diffs$padj <- diff_dev$p_value_adjusted
  list(diffs = diffs,
       scores = scores,
       devObj = dev) -> res
  attr(res,"class") <- "CRCResult"
  res
}


#' plot_cliques
#' @export
plot_cliques <- function(object, plot=c("volcano","tsne","clustergram"),
                         contrast=NULL) {
  stopifnot(class(object) == "CRCResult")
  print('yay')
}

#' @import ggplot2
plot_cliques_volcano <- function(diffs, contrast=NULL) {
  diff_score <- diffs[[contrast[1]]] - diffs[[contrast[2]]]
  toplot <- tibble::tibble(x = diff_score, y= -log10(diffs$padj))
  ggplot(toplot, aes(x=x, y=y, color=x)) + geom_point() +
    theme_classic() + geom_hline(yintercept = -log10(0.05),color='darkgray') +
    scale_color_gradient2(mid = "lightgray", low = "blue",high = "red")
}
