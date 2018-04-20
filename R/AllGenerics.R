# generics --------------------------------------------------------------------
setGeneric("name", function(x) standardGeneric("name"))
setGeneric("name<-", function(x, value) standardGeneric("name<-"))
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
setGeneric("extract_cliques", function(object) standardGeneric("extract_cliques"))
setGeneric("unique_cliques", function(object) standardGeneric("unique_cliques"))
setGeneric("tfbs", function(object) standardGeneric("tfbs"))
setGeneric("hash", function(clique) standardGeneric("hash"))
setGeneric("members", function(clique) standardGeneric("members"))
setGeneric("hashes", function(cliquelist) standardGeneric("hashes"))

# set methods -----------------------------------------------------------------
setMethod("tfbs",signature(object = "CRCView"),
          function(object) {
            object@tfbs
          }
)

setMethod("tfbs",signature(object = "CRCViewList"),
          function(object) {
            object@listData %>% lapply(tfbs) -> nested_tfbs
            all_tf_names <- nested_tfbs %>% lapply(names) %>%
              unlist %>% unique
            top_level_tfbs_by_motif <- list()
            for (i in 1:length(all_tf_names)) {
              all_tf_names[i] -> idx
              nested_tfbs %>% lapply(`[[`,idx) %>% unlist %>%
                GRangesList() %>% unlist %>% reduce -> i_tfbs
              top_level_tfbs_by_motif[[idx]] <- i_tfbs
            }
            GRangesList(top_level_tfbs_by_motif)
          }
)

setMethod("unique_cliques",signature(object = "CliqueList"),
          function(object) {
            names(object) %>% unique %>% object[.]
          }
)


setMethod("extract_cliques",signature(object = "CliqueList"),
          function(object) {
            object@listData
          }
)


setMethod("extract_cliques",signature(object = "CRCView"),
          function(object) {
            object@cliques
          }
)

setMethod("extract_cliques",signature(object = "CRCViewList"),
          function(object) {
            object@listData %>%
              lapply(extract_cliques) %>%
              lapply(extract_cliques) %>%
              CliqueList()
          }
)


setMethod("name",
          signature(x = "TranscriptionFactor"),
          function(x) {
            return(x@name)
          })

setMethod("name",
          signature(x = "CliqueList"),
          function(x) {
            return(x@name)
          })

setMethod("name",
          signature(x = "CRCView"),
          function(x) {
            return(x@name)
          })


setReplaceMethod("name",signature(x="TranscriptionFactor",
                                  value="character"),
                 function(x, value){
                   x@name <- value
                   x
                 }
)

setReplaceMethod("name",signature(x="CliqueList",
                                  value="character"),
                 function(x, value){
                   x@name <- value
                   x
                 }
)


setMethod("pwms",
          signature(x = "TranscriptionFactor"),
          function(x) {
            return(x@pwms)
          })



# one for lists
setMethod("combine",
          signature(x = "list"),
          function(x, ..., name = NULL ) {
            w <- unlist(list(x,...))
            new_name <- name
            if (is.null(name)) new_name <- name(w[[1]])
            pwms <- do.call('c',
                            unlist(lapply(w,
                                          FUN = function(z) pwms(z))))
            return(TranscriptionFactor(name = new_name, pwms = pwms))
          }
)

# one for Regular TranscriptionFactor objects
setMethod("combine",
          signature(x = "TranscriptionFactor"),
          function(x, ..., name = NULL ) {
            w <- unlist(list(x,...))
            new_name <- name
            if (is.null(name)) new_name <- name(w[[1]])
            pwms <- do.call('c',
                            unlist(lapply(w,
                                          FUN = function(z) pwms(z))))
            return(TranscriptionFactor(name = new_name, pwms = pwms))
          }
)


setMethod("hash",
          signature(clique = "Clique"),
          function(clique) {
            return(clique@hash)
          })


setMethod("members",
          signature(clique = "Clique"),
          function(clique) {
            return(clique@members)
          })

setMethod("hashes",
          signature(cliquelist = "CliqueList"),
          function(cliquelist) {
            return(cliquelist@clique_hashes)
          })


 setMethod("show","CRCView", function(object) {
   # print name
   cat("Class: ", class(object),"\n")
   cat("Cliques: ", length(object@cliques), "\n")
   cat("Accessible Sites: ", length(object), "\n")
   cat("TFBS: ", length(head(names(object@tfbs)))," and ",length(object@tfbs)-6," more","\n")
   cat("Bam: ", object@bam, "\n")
 })
