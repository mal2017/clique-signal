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
setGeneric("hash", function(clique) standardGeneric("hash"))
setGeneric("members", function(clique) standardGeneric("members"))
setGeneric("hashes", function(cliquelist) standardGeneric("hashes"))

# set methods -----------------------------------------------------------------
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
   cat("colData: ", !is.null(object@colData), "\n")
   cat("Accessible Sites: ", length(object), "\n")
   cat("TFBS: ", length(object@tfbs), "\n")
   cat("Bam: ", object@bam, "\n")
 })
