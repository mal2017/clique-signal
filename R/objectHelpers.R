#' list_to_cliquelist
#' @export
list_to_cliquelist <- function(cliques) {
  vanilla_clique_to_cliquelist <- function(x) {
    Clique(lapply(x, TranscriptionFactor))
  }
  lapply(cliques, vanilla_clique_to_cliquelist) %>%
    CliqueList()
}
