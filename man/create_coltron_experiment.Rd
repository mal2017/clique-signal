% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importColtron.R
\name{create_coltron_experiment}
\alias{create_coltron_experiment}
\title{Import a set of coltron samples.}
\usage{
create_coltron_experiment(metadata, quantsites = "SUBPEAKS",
  resizeWidth = 1000, quantmode = "READS", nthreads = 1,
  cohort_name = "Coltron", genome = "BSgenome.Hsapiens.UCSC.hg19")
}
\arguments{
\item{metadata}{data.frame with sample metadata.}

\item{quantsites}{Default to subpeaks file in coltron directories. Else paths to bed files.}

\item{resizeWidth}{Width to resize subpeaks to.}

\item{quantmode}{'READS': Full read overlap; 'ATAC': 5' cut sites.}

\item{nthreads}{Number of threads to use for read quant. Passed to Rsubread::featureCounts.}

\item{cohort_name}{Name for experiment group.}

\item{genome}{String name of a BSgenome, used for GC bias addition.}
}
\value{
CRCExperiment object.
}
\description{
User supplies a data.frame containing at least the following
columns:
SAMPLE, CONDITION, COLTRONDIR, BAM.
}
\details{
QUANTSITES is an optional field for specifying atac regions to
search for motifs in and quantify signal. Defaults to coltron's subpeaks
bed file.

QUANTMODE is an optional field for specifying whether to count reads like an atac
or a chipseq experiment. 'ATAC' or 'READS'.

Returns a CRCExperiment, which inherits from RangedSummarizedExperiment and is
compatible with the regular chromVAR workflow as well as clique level analysis.
}
