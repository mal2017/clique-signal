#' Index a CRCExperiment ranges by clique motif presence.
#'
#' @param object A CRCExperiment object.
#' @return A sparse Matrix holding clique motif presence/absence.
#' @export
index_by_clique <- function(object) {
    stopifnot(class(object) == "CRCExperiment")
    cliquewise_tfbs <- tfbs_by_clique(object)
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
