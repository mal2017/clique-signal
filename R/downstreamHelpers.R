#' Index a CRCExperiment ranges by clique motif presence.
#'
#' @param object A CRCExperiment object.
#' @param remove_subsets Bool. Remove cliques that are subsets of other cliques.
#' @param combine_similar FALSE or int n. Recursively combine cliques until all cliques have more than n differences from all others.
#' @return A sparse Matrix holding clique motif presence/absence.
#' @export
index_by_clique <- function(object, remove_subsets = T, combine_similar = F) {
    stopifnot(class(object) == "CRCExperiment")
    cliquewise_tfbs <- tfbs_by_clique(object,
                                      remove_subsets = remove_subsets,
                                      combine_similar = combine_similar)
    # https://support.bioconductor.org/p/82452/
    message("Indexing accessible sites by clique-grouped TFBS...", appendLF = F)
    tictoc::tic()
    hits <- findOverlaps(rowRanges(object), cliquewise_tfbs)
    mat <- matrix(FALSE, nrow = S4Vectors::queryLength(hits), ncol = S4Vectors::subjectLength(hits))
    mat[S4Vectors::as.matrix(hits)] <- TRUE
    colnames(mat) <- names(cliquewise_tfbs)
    rownames(mat) <- rownames(object)
    tictoc::toc()
    Matrix::Matrix(mat, sparse = T)
}

#' Run differential analysis at the clique level.
#'
#' Supply a CRCExperiment, and index matrix.
#'
#' See \link[chromVAR]{computeDeviations}.
#' See \link[chromVAR]{differentialDeviations}.
#'
#' @param object A CRCExperiment object.
#' @param ix A logical matrix holding clique index information.
#' @param contrast Optionally provide vector containing 2 levels of metadata CONDITION.
#' @param parametric Passed to \link[chromVAR]{differentialDeviations}.
#' @param score_type Use chromVAR's deviations or deviation z-scores as metric for results table.
#' @return S3 object of class CRCResult, a list holding means+pvalues, sample scores, and a chromVARdeviations object.
#' @export
get_diff_cliques <- function(object, ix, contrast = NULL, parametric = TRUE, score_type = c("dev",
    "z")) {
    score_type <- score_type[1]
    stopifnot(nrow(object) == nrow(ix))
    stopifnot(length(contrast) == 2 | is.null(contrast))

    if (!is.null(contrast)) {
        object <- object[, object$CONDITION %in% contrast]
    }

    colData(object) <- droplevels(colData(object))
    message("Computing deviations via chromVAR... ", appendLF = F)
    tictoc::tic()
    dev <- chromVAR::computeDeviations(object, annotations = ix)
    diff_dev <- chromVAR::differentialDeviations(dev, groups = "CONDITION", alternative = "two.sided",
        parametric = parametric)
    tictoc::toc()

    if (score_type == "dev") {
        scores <- chromVAR::deviations(dev)
        scores <- tibble::as_tibble(scores, rownames = "clique")
    } else if (score_type == "z") {
        scores <- chromVAR::deviationScores(dev)
        scores <- tibble::as_tibble(scores, rownames = "clique")
    }

    names_by_cond <- split(rownames(colData(object)), colData(object)$CONDITION)

    diffs <- tibble::tibble(clique = scores$clique)
    for (i in 1:length(names_by_cond)) {
        grp_name <- names(names_by_cond)[i]
        samples <- names_by_cond[[i]]
        diffs[[grp_name]] <- rowMeans(scores[samples])
    }

    diffs$padj <- diff_dev$p_value_adjusted
    diffs$pval <- diff_dev$p_value
    res <- list(diffs = diffs, scores = scores, devObj = dev)
    attr(res, "class") <- "CRCResult"
    res
}


#' Plotting helper for CRCResult object.
#'
#' @param object An object of class CRCResult.
#' @param plot A string specifying the type of plot to produce. Only 'volcano' implemented for now.
#' @param contrast A vector containing the 2 levels of CONDITION to compare.
#' @param use.adjusted.p Logical specifying whether to use adjusted or non-adjusted p-values for plotting.
#' @return ggplot2 object.
#' @export
plot_cliques <- function(object, plot = c("volcano", "tsne", "clustergram"), contrast = NULL, use.adjusted.p = F) {
    stopifnot(class(object) == "CRCResult")
    if (plot == "volcano") {
        if (use.adjusted.p) {
            p <- "padj"
        } else {
            p <- "pval"
        }
        plot_cliques_volcano(object$diffs, contrast = contrast, p = p)
    } else {
        message("Only volcano helper available in this version.")
    }
}

#' @import ggplot2
plot_cliques_volcano <- function(diffs, p, contrast = NULL) {
    diff_score <- diffs[[contrast[1]]] - diffs[[contrast[2]]]
    toplot <- tibble::tibble(x = diff_score, y = -log10(diffs[[p]]))
    lims_x <- abs(toplot$x) %>% max %>% c(-1 * ., .)
    lims_y <- c(0, 1.2 * max(toplot$y))
    ggplot(toplot, aes(x = x, y = y, color = x)) + geom_point() + theme_classic() + geom_hline(yintercept = -log10(0.05),
        color = "darkgray") + scale_color_gradient2(mid = "lightgray", low = "blue", high = "red") +
        xlim(lims_x) + ylim(lims_y)
}


#' Remove vecs that are perfect subsets of other vecs
#'
#' Not really intended to be called by user.
#'
#' @param named_list_of_vectors list of vecs
#' @import foreach
#' @import doParallel
remove_subset_vectors <- function(named_list_of_vectors) {
  obj <- unique(named_list_of_vectors)
  names_obj <- lapply(obj, paste, collapse = ",")
  names(obj) <- names_obj
  seq <- 1:length(obj)
  ix <- foreach(i = seq, .combine = 'c') %dopar% {
    subsets_of_i <- c()
    for (k in seq) {
      if (i == k) next
      if (all(obj[[k]] %in% obj[[i]])) {
        subsets_of_i <- c(subsets_of_i, k)
      }
    }
    subsets_of_i
  }

  obj[-unique(ix)]
}

# ------------------------------

get_merged_vec <- function(vec, list_of_vecs, max_diff = NULL) {
  diffs <- lapply(list_of_vecs, get_setdiff_total, vec)
  which(diffs <= max_diff) -> ix
  list_of_vecs[ix] %>% unlist %>% unique %>% sort -> new_vec
  new_vec
}

get_setdiff_total <- function(a,b) {
  length(setdiff(a,b)) + length(setdiff(b,a))
}

#' Combine similar vecs
#'
#' Not really intended to be called by user.
#'
#' Uses doparallel backend, register this with
#' 'registerDoParallel(12)'.
#'
#' @param named_list_of_vectors list of vecs
#' @param max_diff_to_combine recursively combine until all cliques are have more diffs than this arg.
#' @param verbose print progress for debugging
#' @import foreach
#' @import doParallel
combine_similar_vectors <- function(named_list_of_vectors,
                                     max_diff_to_combine = 1, verbose =F) {
  new_obj <- foreach(v = named_list_of_vectors) %dopar% {
    get_merged_vec(v, named_list_of_vectors,
                   max_diff = max_diff_to_combine)
  }
  names_new_obj <- lapply(new_obj, paste, collapse = ",") # sorted above
  names(new_obj) <- names_new_obj
  new_obj %<>%
    names %>%
    sort %>%
    duplicated %>% not %>%
    which %>% new_obj[.]
  return(new_obj)
}
