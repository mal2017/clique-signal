# generics --------------------------------------------------------------------

#' @rdname name
#' @export
setGeneric("name", function(x) standardGeneric("name"))

#' @rdname name
#' @export
setGeneric("name<-", function(x, value) standardGeneric("name<-"))

#' @rdname pwms
#' @export
setGeneric("pwms", function(x) standardGeneric("pwms"))
setGeneric("combine",
           function(x, ..., name=name(x)) {
             if (length(list(...)) > 1) {
               combine(x, do.call(combine, list(...)))
             } else {
               standardGeneric("combine")
             }
           }
)

#' @rdname extract_cliques
#' @export
setGeneric("extract_cliques",
           function(object) standardGeneric("extract_cliques"))

#' @rdname unique_cliques
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

#' @rdname members
#' @export
setGeneric("members", function(clique) standardGeneric("members"))

#' @rdname name
#' @export
setGeneric("hashes", function(cliquelist) standardGeneric("hashes"))

#' @rdname bam
#' @export
setGeneric("bam", function(object) standardGeneric("bam"))

# set methods -----------------------------------------------------------------

#' @name bam
#' @rdname CRCViewList
#' @export
setMethod("bam", signature(object = "CRCViewList"),
          function(object) {
            object %>% lapply(bam)
          }
)

#' @name bam
#' @rdname CRCView
#' @export
setMethod("bam", signature(object = "CRCView"),
          function(object) {
            object@bam
          }
)

#' @name quantsites
#' @rdname CRCView
#' @export
setMethod("quantsites", signature(object = "CRCView"),
          function(object) {
            GRanges(object)
          }
)

#' @name quantsites
#' @rdname CRCViewList
#' @export
setMethod("quantsites", signature(object = "CRCViewList"),
          function(object) {
            object@listData %>% lapply(GRanges) %>% GRangesList
          }
)

#' @name tfbs
#' @rdname CRCView
#' @export
setMethod("tfbs", signature(object = "CRCView"),
          function(object) {
            object@tfbs
          }
)

#' @name tfbs
#' @rdname CRCViewList
#' @export
setMethod("tfbs", signature(object = "CRCViewList"),
          function(object) {
            object@listData %>% lapply(tfbs) -> nested_tfbs
            all_tf_names <- nested_tfbs %>% lapply(names) %>%
              unlist %>% unique
            top_level_tfbs_by_motif <- list()
            for (i in 1:length(all_tf_names)) {
              all_tf_names[i] -> idx
              nested_tfbs %>% lapply(`[[`, idx) %>% unlist %>%
                GRangesList() %>% unlist %>%
                GenomicRanges::reduce() -> i_tfbs
              top_level_tfbs_by_motif[[idx]] <- i_tfbs
            }
            GRangesList(top_level_tfbs_by_motif)
          }
)

#' @name tfbs
#' @rdname CRCExperiment
#' @export
setMethod("tfbs", signature(object = "CRCExperiment"),
          function(object) {
            tfbs(object@crcs)
          }
)

#' @name unique_cliques
#' @rdname CliqueList
#' @export
setMethod("unique_cliques", signature(object = "CliqueList"),
          function(object) {
            names(object) %>% unique %>% object[.]
          }
)

#' @name extract_cliques
#' @rdname CliqueList
#' @export
setMethod("extract_cliques", signature(object = "CliqueList"),
          function(object) {
            object@listData
          }
)

#' @name extract_cliques
#' @rdname CRCView
#' @export
setMethod("extract_cliques", signature(object = "CRCView"),
          function(object) {
            object@cliques
          }
)

#' @name extract_cliques
#' @rdname CRCViewList
#' @export
setMethod("extract_cliques", signature(object = "CRCViewList"),
          function(object) {
            object@listData %>%
              lapply(extract_cliques) %>%
              lapply(extract_cliques) %>%
              CliqueList()
          }
)

#' @name extract_cliques
#' @rdname CRCExperiment
#' @export
setMethod("extract_cliques", signature(object = "CRCExperiment"),
          function(object) {
            extract_cliques(object@crcs)
          }
)

#' @name name
#' @rdname TranscriptionFactor
#' @export
setMethod("name",
          signature(x = "TranscriptionFactor"),
          function(x) {
            return(x@name)
          })

#' @name name
#' @rdname CliqueList
#' @export
setMethod("name",
          signature(x = "CliqueList"),
          function(x) {
            return(x@name)
          })

#' @name name
#' @rdname CRCView
#' @export
setMethod("name",
          signature(x = "CRCView"),
          function(x) {
            return(x@name)
          })

#' @name name
#' @rdname TranscriptionFactor
#' @exportMethod "name<-"
setReplaceMethod("name", signature(x = "TranscriptionFactor",
                                  value = "character"),
                 function(x, value){
                   x@name <- value
                   x
                 }
)

#' @name name
#' @rdname CliqueList
#' @exportMethod "name<-"
setReplaceMethod("name", signature(x = "CliqueList",
                                  value = "character"),
                 function(x, value){
                   x@name <- value
                   x
                 }
)

#' @name pwms
#' @rdname TranscriptionFactor
#' @export
setMethod("pwms",
          signature(x = "TranscriptionFactor"),
          function(x) {
            return(x@pwms)
          })



#' @name combine
#' @rdname TranscriptionFactor
#' @export
setMethod("combine",
          signature(x = "list"),
          function(x, ..., name = NULL ) {
            w <- unlist(list(x, ...))
            new_name <- name
            if (is.null(name)) new_name <- name(w[[1]])
            pwms <- do.call("c",
                            unlist(lapply(w,
                                          FUN = function(z) pwms(z))))
            return(TranscriptionFactor(name = new_name, pwms = pwms))
          }
)

#' @name combine
#' @rdname TranscriptionFactor
#' @export
setMethod("combine",
          signature(x = "TranscriptionFactor"),
          function(x, ..., name = NULL ) {
            w <- unlist(list(x, ...))
            new_name <- name
            if (is.null(name)) new_name <- name(w[[1]])
            pwms <- do.call("c",
                            unlist(lapply(w,
                                          FUN = function(z) pwms(z))))
            return(TranscriptionFactor(name = new_name, pwms = pwms))
          }
)


#' @name hash
#' @rdname Clique
#' @export
setMethod("hash",
          signature(clique = "Clique"),
          function(clique) {
            return(clique@hash)
          })

#' @name members
#' @rdname Clique
#' @export
setMethod("members",
          signature(clique = "Clique"),
          function(clique) {
            return(clique@members)
          })

#' @name hashes
#' @rdname CliqueList
#' @export
setMethod("hashes",
          signature(cliquelist = "CliqueList"),
          function(cliquelist) {
            return(cliquelist@clique_hashes)
          })


#' @name show
#' @rdname CRCView
#' @export
setMethod("show", "CRCView", function(object) {
  # print name
  cat("Class: ", class(object), "\n")
  cat("Cliques: ", length(object@cliques), "\n")
  cat("Accessible Sites: ", length(object), "\n")
  cat("TFBS: ", length(head(names(object@tfbs))), " and ",
      length(object@tfbs) - 6, " more", "\n")
  cat("Bam: ", object@bam, "\n")
})

#' @name show
#' @rdname CRCExperiment
#' @export
setMethod("show", "CRCExperiment", function(object) {
  # print name
  cat("Class: ", class(object), "\n")
  cat("Samples: ", length(object@crcs), "\n")
  cat("Accessible Sites: ", length(object), "\n")
  cat("Conditions: ", crc@metadata$CONDITION %>% levels, "\n")
})
