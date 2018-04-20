#' list_to_cliquelist
#' @export
list_to_cliquelist <- function(cliques) {
  vanilla_clique_to_cliquelist <- function(x) {
    Clique(lapply(x, TranscriptionFactor))
  }
  lapply(cliques, vanilla_clique_to_cliquelist) %>%
    CliqueList()
}

#' tfbs_by_clique
#' @export
tfbs_by_clique <- function(object) {
  tfbs_by_tf <- tfbs(object)
  cliques <- unique_cliques(extract_cliques(object))
  get_tfbs_gr <- function(tfnames,gr=tfbs_by_tf) {
    gr[tfnames] %>% unlist %>% reduce
  }
  message("Grouping TFBS by Clique...", appendLF = F)
  tictoc::tic()
  lapply(cliques,members) %>%
    lapply(get_tfbs_gr) %>%
    GRangesList() -> res
  tictoc::toc()
  res
}
