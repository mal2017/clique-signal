% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downstreamHelpers.R
\name{index_by_clique}
\alias{index_by_clique}
\title{Index a CRCExperiment ranges by clique motif presence.}
\usage{
index_by_clique(object, remove_subsets = T, combine_similar = F)
}
\arguments{
\item{object}{A CRCExperiment object.}

\item{remove_subsets}{Bool. Remove cliques that are subsets of other cliques.}

\item{combine_similar}{FALSE or int n. Recursively combine cliques until all cliques have more than n differences from all others.}
}
\value{
A sparse Matrix holding clique motif presence/absence.
}
\description{
Index a CRCExperiment ranges by clique motif presence.
}
