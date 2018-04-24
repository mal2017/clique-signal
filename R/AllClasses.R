
setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("GRangesListOrNULL", c("GRangesList", "NULL"))
setClassUnion("GRangesOrNULL", c("GRanges", "NULL"))
setClassUnion("DataFrameOrNULL", c("DataFrame", "NULL"))

#' A class for holding Transcription Factors.
#'
#' Holds a name and a PWMatrixList, which can be empty.
#'
#' @rdname TranscriptionFactor
#' @slot name Name for the TF
#' @slot pwms PWMatrixList. Defined in TFBSTools.
#' @export
setClass("TranscriptionFactor", representation(name = "character", pwms = "PWMatrixList"))

#' @rdname TranscriptionFactor
#'
#' @param name A name for the TranscriptionFactor.
#' @param pwms Optionally a PWMatrixList object.
#' @export
TranscriptionFactor <- function(name, pwms = NULL) {
    if (is.null(pwms)) {
        pwmlist <- TFBSTools::PWMatrixList()
    } else if (class(pwms) == "PWMatrixList") {
        pwmlist <- pwms
    } else {
        pwmlist <- do.call("c", lapply(c(pwms), TFBSTools::PWMatrixList))
    }
    new("TranscriptionFactor", name = name, pwms = pwmlist)
}

# -----------------------------------------------------------------------------

#' A class for holding groups of TranscriptionFactor objects.
#'
#' Contains at least one TranscriptionFactor object, as well as
#' a list of names of contained TF instances.
#'
#' @rdname Clique
#' @export
setClass("Clique", contains = "TranscriptionFactor", representation(members = "character", hash = "ANY"))

#' @rdname Clique
#'
#' @param tf TranscriptionFactor object.
#' @param ... more TranscriptionFactor objects.
#' @param name A name for the clique.
#' @export
Clique <- function(tf, ..., name = NULL) {
    others <- list(...)
    tfs <- c(tf, others)
    members <- sort(unlist(lapply(tfs, FUN = function(x) name(x))))
    derived_name <- paste(members, collapse = ",")
    hash <- substr(openssl::md5(derived_name), 1, 8)
    new_clique <- do.call(combine, tfs)
    if (is.null(name)) {
        name(new_clique) <- derived_name
    } else {
        name(new_clique) <- name
    }
    new("Clique", new_clique, members = members, hash = hash)
}

# -----------------------------------------------------------------------------

#' A class for holding groups of Clique objects.
#'
#' @rdname CliqueList
#' @export
setClass("CliqueList", contains = "SimpleList", representation(name = "characterOrNULL", clique_names = "character",
    clique_hashes = "ANY"))

#' @rdname CliqueList
#' @param clique Clique object.
#' @param ... more Clique objects.
#' @param name Sample name.
#' @export
CliqueList <- function(clique, ..., name = NULL) {
    cliques <- unlist(list(clique, ...))
    clique_names <- unlist(lapply(cliques, name))
    clique_hashes <- unlist(lapply(cliques, hash))
    new_clqs <- new("CliqueList", SimpleList(cliques), name = name, clique_names = clique_names,
        clique_hashes = clique_hashes)
    names(new_clqs) <- clique_names
    new_clqs
}

# -----------------------------------------------------------------------------

#' A class for holding a representation of a single sample's predicted
#' network.
#'
#' Contains a GRanges object holding accessible regions for downstream consideration.
#' Additionally may contain a name, a CliqueList, all predicted TFBS, and a path
#' to a bam used for quantification.
#'
#' @rdname CRCView
#' @export
setClass("CRCView", contains = "GRanges", representation(name = "characterOrNULL", cliques = "CliqueList",
    tfbs = "GRangesListOrNULL", bam = "characterOrNULL"))

#' @rdname CRCView
#' @param quantsites GRanges object.
#' @param cliques CliqueList object.
#' @param prior_tfbs Optional GrangesList of predicted TFBS.
#' @param sample Sample name.
#' @param bampath Path to bam for quantification.
#' @export
CRCView <- function(quantsites, cliques, prior_tfbs = NULL, sample = name(cliques), bampath = NULL) {
    new("CRCView", quantsites, cliques = cliques, name = sample, tfbs = prior_tfbs, bam = bampath)
}

# -----------------------------------------------------------------------------

#' A class for holding a list of CRCViews.
#'
#' @rdname CRCViewList
#' @export
setClass("CRCViewList", contains = "SimpleList", representation(samples = "characterOrNULL"))

#' @rdname CRCViewList
#' @param crcview A CRCView object.
#' @param ... More CRCView Objects.
#' @param name A name for the cohort.
#' @export
CRCViewList <- function(crcview, ..., name = NULL) {
    crcs <- unlist(list(crcview, ...))
    crc_names <- unlist(lapply(crcs, name))
    new_crcv <- new("CRCViewList", SimpleList(crcs), samples = crc_names)
    names(new_crcv) <- crc_names
    new_crcv
}

# -----------------------------------------------------------------------------

#' A class for containing the full representation of a cohort's coltron
#' results.
#'
#' @rdname CRCExperiment
#' @export
setClass("CRCExperiment", contains = "RangedSummarizedExperiment", representation(name = "characterOrNULL",
    crcs = "CRCViewList"))

#' @rdname CRCExperiment
#'
#' @param rse A RangedSummarizedExperiment object.
#' @param crclist A CRCList object.
#' @param cohort Optional name for experimental group.
#' @export
CRCExperiment <- function(rse, crclist, cohort = NULL) {
    new("CRCExperiment", rse, crcs = crclist, name = cohort)
}
