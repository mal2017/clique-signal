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
#' @param combine_similar FALSE or int n. Recursively combine cliques until all cliques have more than n differences from all others.
#' @export
tfbs_by_clique <- function(object, remove_subsets = T, combine_similar = F) {
    tfbs_by_tf <- tfbs(object)
    cliques <- unique_cliques(extract_cliques(object))
    if (remove_subsets) cliques <- remove_subset_cliques(cliques)
    if (combine_similar) cliques <- combine_similar_cliques(cliques,
                                                            combine_when_at_least_as_similar = combine_similar)
    get_tfbs_gr <- function(tfnames, gr = tfbs_by_tf) {
        gr[tfnames] %>% unlist %>% GenomicRanges::reduce()
    }
    message("Grouping TFBS by Clique... ", appendLF = F)
    tictoc::tic()
    res <- lapply(cliques, members) %>% lapply(get_tfbs_gr) %>% GRangesList()
    tictoc::toc()
    res
}
