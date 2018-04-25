# -----------------------------------------------------------------------------
# generics --------------------------------------------------------------------
# -----------------------------------------------------------------------------

#' @rdname name
#' @export
setGeneric("name", function(x) standardGeneric("name"))

#' @rdname name
#' @export
setGeneric("name<-", function(x, value) standardGeneric("name<-"))

#' @rdname pwms
#' @export
setGeneric("pwms", function(x) standardGeneric("pwms"))

#' @rdname combine
#' @export
setGeneric("combine", function(x, ..., name = name(x)) {
    if (length(list(...)) > 1) {
        combine(x, do.call(combine, list(...)))
    } else {
        standardGeneric("combine")
    }
})

#' @rdname extract_cliques
#' @export
setGeneric("extract_cliques", function(object) standardGeneric("extract_cliques"))

#' @rdname extract_cliques
#' @export
setGeneric("unique_cliques", function(object) standardGeneric("unique_cliques"))

#' @rdname tfbs
#' @export
setGeneric("tfbs", function(object) standardGeneric("tfbs"))

#' @rdname hash
#' @export
setGeneric("hash", function(clique) standardGeneric("hash"))

#' @rdname quantsites
#' @export
setGeneric("quantsites", function(object) standardGeneric("quantsites"))

#' @rdname name
#' @export
setGeneric("members", function(clique) standardGeneric("members"))

#' @rdname hash
#' @export
setGeneric("hashes", function(cliquelist) standardGeneric("hashes"))

#' @rdname bam
#' @export
setGeneric("bam", function(object) standardGeneric("bam"))

#' @rdname show
#' @export
setGeneric("show", function(object) standardGeneric("show"))

#' @rdname remove_subset_cliques
#' @export
setGeneric("remove_subset_cliques", function(object) standardGeneric("remove_subset_cliques"))

# set methods -----------------------------------------------------------------

#' Remove cliques that are subsets of other cliques
#' @rdname remove_subset_cliques
#' @export
setMethod("remove_subset_cliques", signature(object = "CliqueList"),
          function(object) {
            not_subsets <- lapply(object, members) %>%
              remove_subset_vectors() %>% names
            object[not_subsets]
})

#' @rdname remove_subset_cliques
#' @export
setMethod("remove_subset_cliques", signature(object = "CRCExperiment"),
          function(object) {
            not_subsets <- lapply(object, members) %>%
              remove_subset_vectors() %>% names
            object[not_subsets]
          })


# -----------------------------------------------------------------------------

#' Get bam
#' @rdname bam
#' @param object Object holding bam paths.
#' @export
setMethod("bam", signature(object = "CRCViewList"), function(object) {
    object %>% lapply(bam)
})

#' @rdname bam
#' @export
setMethod("bam", signature(object = "CRCView"), function(object) {
    object@bam
})

# -----------------------------------------------------------------------------

#' Get quantsites.
#' @rdname quantsites
#' @param object Object holding quantsites.
#' @export
setMethod("quantsites", signature(object = "CRCView"), function(object) {
    GRanges(object)
})

#' @rdname quantsites
#' @export
setMethod("quantsites", signature(object = "CRCViewList"), function(object) {
    object@listData %>% lapply(GRanges) %>% GRangesList
})

# -----------------------------------------------------------------------------

#' Get tfbs.
#' @rdname tfbs
#' @param object Object holding tfbs as a GRanges object.
#' @export
setMethod("tfbs", signature(object = "CRCView"), function(object) {
    object@tfbs
})

#' @rdname tfbs
#' @export
setMethod("tfbs", signature(object = "CRCViewList"), function(object) {
    nested_tfbs <- object@listData %>% lapply(tfbs)
    all_tf_names <- nested_tfbs %>% lapply(names) %>% unlist %>% unique
    top_level_tfbs_by_motif <- list()
    for (i in 1:length(all_tf_names)) {
        idx <- all_tf_names[i]
        i_tfbs <- nested_tfbs %>% lapply(`[[`, idx) %>% unlist %>% GRangesList() %>% unlist %>% GenomicRanges::reduce()
        top_level_tfbs_by_motif[[idx]] <- i_tfbs
    }
    GRangesList(top_level_tfbs_by_motif)
})

#' @rdname tfbs
#' @export
setMethod("tfbs", signature(object = "CRCExperiment"), function(object) {
    tfbs(object@crcs)
})

# -----------------------------------------------------------------------------

#' Clique extractors.
#' @rdname extract_cliques
#' @param object Object holding cliques.
#' @export
setMethod("unique_cliques", signature(object = "CliqueList"), function(object) {
    names(object) %>% unique %>% object[.]
})

#' @rdname extract_cliques
#' @export
setMethod("extract_cliques", signature(object = "CliqueList"), function(object) {
    object@listData
})

#' @rdname extract_cliques
#' @export
setMethod("extract_cliques", signature(object = "CRCView"), function(object) {
    object@cliques
})

#' @rdname extract_cliques
#' @export
setMethod("extract_cliques", signature(object = "CRCViewList"), function(object) {
    object@listData %>% lapply(extract_cliques) %>% lapply(extract_cliques) %>% CliqueList()
})

#' @rdname extract_cliques
#' @export
setMethod("extract_cliques", signature(object = "CRCExperiment"), function(object) {
    extract_cliques(object@crcs)
})

# -----------------------------------------------------------------------------

#' Name getters.
#' @rdname name
#' @param x A name-able object.
#' @param value New name.
#' @export
setMethod("name", signature(x = "TranscriptionFactor"), function(x) {
    return(x@name)
})

#' @rdname name
#' @export
setMethod("name", signature(x = "CliqueList"), function(x) {
    return(x@name)
})

#' @rdname name
#' @export
setMethod("name", signature(x = "CRCView"), function(x) {
    return(x@name)
})

#' @rdname name
#' @exportMethod 'name<-'
setReplaceMethod("name", signature(x = "TranscriptionFactor", value = "character"), function(x, value) {
    x@name <- value
    x
})

#' @rdname name
#' @exportMethod 'name<-'
setReplaceMethod("name", signature(x = "CliqueList", value = "character"), function(x, value) {
    x@name <- value
    x
})

#' @rdname name
#' @param clique A clique object.
#' @export
setMethod("members", signature(clique = "Clique"), function(clique) {
  return(clique@members)
})

# -----------------------------------------------------------------------------

#' Get pwms.
#' @rdname pwms
#' @param x A TranscriptionFactor object.
#' @export
setMethod("pwms", signature(x = "TranscriptionFactor"), function(x) {
    return(x@pwms)
})

# -----------------------------------------------------------------------------

#' Combine TranscriptionFactor instances.
#' @rdname combine
#' @param x A TranscriptionFactor object or list of TranscriptionFactor objects.
#' @param ... TranscriptionFactor objects.
#' @param name A new name for the combined TF. Optional. Defaults to name of first TF.
#' @export
setMethod("combine", signature(x = "list"), function(x, ..., name = NULL) {
    w <- unlist(list(x, ...))
    new_name <- name
    if (is.null(name))
        new_name <- name(w[[1]])
    pwms <- do.call("c", unlist(lapply(w, FUN = function(z) pwms(z))))
    return(TranscriptionFactor(name = new_name, pwms = pwms))
})

#' @rdname combine
#' @export
setMethod("combine", signature(x = "TranscriptionFactor"), function(x, ..., name = NULL) {
    w <- unlist(list(x, ...))
    new_name <- name
    if (is.null(name))
        new_name <- name(w[[1]])
    pwms <- do.call("c", unlist(lapply(w, FUN = function(z) pwms(z))))
    return(TranscriptionFactor(name = new_name, pwms = pwms))
})

# -----------------------------------------------------------------------------

#' Get unique hash.
#' @rdname hash
#' @param clique A Clique.
#' @param cliquelist A CliqueList.
#' @export
setMethod("hash", signature(clique = "Clique"), function(clique) {
    return(clique@hash)
})

#' @rdname hash
#' @export
setMethod("hashes", signature(cliquelist = "CliqueList"), function(cliquelist) {
  return(cliquelist@clique_hashes)
})

# -----------------------------------------------------------------------------

#' Print object representation.
#' @rdname show
#' @param object Object with printable representation.
#' @export
setMethod("show", "CRCView", function(object) {
    # print name
    cat("Class: ", class(object), "\n")
    cat("Cliques: ", length(object@cliques), "\n")
    cat("Accessible Sites: ", length(object), "\n")
    cat("TFBS: ", length(head(names(object@tfbs))), " and ", length(object@tfbs) - 6, " more", "\n")
    cat("Bam: ", object@bam, "\n")
})

#' @rdname show
#' @export
setMethod("show", "CRCExperiment", function(object) {
    # print name
    cat("Class: ", class(object), "\n")
    cat("Samples: ", length(object@crcs), "\n")
    cat("Accessible Sites: ", length(object), "\n")
    cat("Conditions: ", object@metadata$CONDITION %>% levels, "\n")
})
