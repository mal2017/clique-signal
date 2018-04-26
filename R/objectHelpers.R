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
#' @param combine_similar FALSE or int n. Combine cliques with less than N differences.
#' @import foreach
#' @import doParallel
#' @export
tfbs_by_clique <- function(object, remove_subsets = T, combine_similar = F) {
  message("Extracting all TFBS and cliques from your experiment... ", appendLF = F)
  tictoc::tic()
  tfbs_by_tf <- tfbs(object)
  cliques <- unique_cliques(extract_cliques(object))
  tictoc::toc()
  if (remove_subsets) cliques <- remove_subset_cliques(cliques)
  if (combine_similar) cliques <- combine_similar_cliques(cliques,
                                                          combine_when_at_least_as_similar = combine_similar)
  get_tfbs_gr <- function(tfnames, gr = tfbs_by_tf) {
    gr[tfnames] %>% unlist %>% GenomicRanges::reduce()
  }
  message("Grouping TFBS by Clique... ", appendLF = F)
  tictoc::tic()
  #res <- lapply(cliques, members) %>% lapply(get_tfbs_gr) %>% GRangesList()
  res <- foreach(x = lapply(cliques, members), .final = function(z) GRangesList(z)) %dopar% {
    get_tfbs_gr(x)
  }
  names(res) <- names(cliques)
  tictoc::toc()
  res
}
