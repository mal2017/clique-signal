# unionClasses ----------------------------------------------------------------
setClassUnion("characterOrNULL",c("character","NULL"))
setClassUnion("GRangesListOrNULL",c("GRangesList","NULL"))
setClassUnion("GRangesOrNULL",c("GRanges","NULL"))
setClassUnion("DataFrameOrNULL",c("DataFrame","NULL"))

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# classes ---------------------------------------------------------------------
# -----------------------------------------------------------------------------

#' @rdname TranscriptionFactor
#' @export
setClass("TranscriptionFactor",
         representation(name = "character",
                        pwms = "PWMatrixList"))

TranscriptionFactor <- function(name, pwms = NULL) {
  if (is.null(pwms)) {
    pwmlist <- TFBSTools::PWMatrixList()
  } else if (class(pwms) == "PWMatrixList") {
    pwmlist <- pwms
  } else {
    pwmlist <- do.call('c',lapply(c(pwms),TFBSTools::PWMatrixList))
  }
  new("TranscriptionFactor",name = name, pwms = pwmlist)
}







# -----------------------------------------------------------------------------

#' @rdname Clique
#' @export
setClass("Clique",
         contains = "TranscriptionFactor",
         representation(members = "character",
                        hash = "ANY"))

Clique <- function(tf, ..., name = NULL) {
  others <- list(...)
  tfs <- c(tf,others)
  members <- sort(unlist(lapply(tfs, FUN = function(x) name(x))))
  derived_name <- paste(members,collapse = ",")
  hash <- substr(openssl::md5(derived_name),1,8)
  new_clique <- do.call(combine,tfs)
  if (is.null(name)) {
    name(new_clique) <- derived_name
  } else {
    name(new_clique) <- name
  }
  new("Clique", new_clique, members = members, hash = hash)
}


# -----------------------------------------------------------------------------

#' @rdname CliqueList
#' @export
setClass("CliqueList",
         contains = "SimpleList",
         representation(name = "characterOrNULL",
                       clique_names = "character",
                        clique_hashes = "ANY")
)

CliqueList <- function(clique, ..., name=NULL ) {
  cliques <- unlist(list(clique, ...))
  clique_names <- unlist(lapply(cliques,name))
  clique_hashes <- unlist(lapply(cliques,hash))
  new("CliqueList", SimpleList(cliques),
      name=name, clique_names = clique_names,
      clique_hashes = clique_hashes) -> new_clqs
  names(new_clqs) <- clique_names
  new_clqs
}



# -----------------------------------------------------------------------------

##' @rdname CRCView
##' @export
setClass("CRCView",
         contains = "GRanges",
         representation(name = "characterOrNULL",
                        cliques = "CliqueList",
                        tfbs = "GRangesListOrNULL",
                        bam = "characterOrNULL"))

CRCView <- function(quantsites, cliques, prior_tfbs=NULL,
                    sample=name(cliques), bampath=NULL){
  new("CRCView", quantsites, cliques=cliques,
      name=sample, tfbs = prior_tfbs, bam=bampath)
}


# -----------------------------------------------------------------------------

#' @rdname CRCViewList
#' @export
setClass("CRCViewList",
         contains = "SimpleList",
         representation(samples = "characterOrNULL"))

CRCViewList <- function(crcview, ..., name=NULL ) {
  crcs <- unlist(list(crcview, ...))
  crc_names <- unlist(lapply(crcs,name))
  new("CRCViewList", SimpleList(crcs),
      samples = crc_names) -> new_crcv
  names(new_crcv) <- crc_names
  new_crcv
}

# -----------------------------------------------------------------------------

#' #' @rdname CRCExperiment
#' #' @export
 setClass("CRCExperiment",
          contains = "RangedSummarizedExperiment",
          representation(name = "characterOrNULL",
                         crcs = "CRCViewList"))

CRCExperiment <- function(rse, crclist, cohort = NULL) {
   new("CRCExperiment", rse, crcs = crclist, name=cohort)
}
