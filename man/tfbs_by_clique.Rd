% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objectHelpers.R
\name{tfbs_by_clique}
\alias{tfbs_by_clique}
\title{Group TFBS  by clique.}
\usage{
tfbs_by_clique(object, remove_subsets = T, combine_similar = F)
}
\arguments{
\item{object}{A CRCExperiment, CRCViewList, or CRCView object.}

\item{remove_subsets}{Bool. Exclude cliques that are subsets of other cliques.}

\item{combine_similar}{FALSE or int n. Combine cliques with less than N differences.}
}
\description{
Group TFBS  by clique.
}
