#' Convert list of lists containing TF names to a CliqueList.
#'
#' @param cliques A list of character vectors.
#' @export
list_to_cliquelist <- function(cliques) {
    vanilla_clique_to_cliquelist <- function(x) {
        Clique(lapply(x, TranscriptionFactor))
    }
    lapply(cliques, vanilla_clique_to_cliquelist) %>% CliqueList()
}

#' Group TFBS  by clique.
#'
#' @param object A CRCExperiment, CRCViewList, or CRCView object.
#' @param remove_subsets Bool. Exclude cliques that are subsets of other cliques.
#' @export
tfbs_by_clique <- function(object, remove_subsets = T) {
    tfbs_by_tf <- tfbs(object)
    cliques <- unique_cliques(extract_cliques(object))
    if (remove_subsets) cliques <- remove_subset_cliques(cliques)
    get_tfbs_gr <- function(tfnames, gr = tfbs_by_tf) {
        gr[tfnames] %>% unlist %>% GenomicRanges::reduce()
    }
    message("Grouping TFBS by Clique... ", appendLF = F)
    tictoc::tic()
    res <- lapply(cliques, members) %>% lapply(get_tfbs_gr) %>% GRangesList()
    tictoc::toc()
    res
}
