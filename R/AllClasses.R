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

setGeneric("name", function(tf) standardGeneric("name"))
setMethod("name",
          signature(tf = "TranscriptionFactor"),
          function(tf) {
            return(tf@name)
          })

setGeneric("pwms", function(tf) standardGeneric("pwms"))
setMethod("pwms",
          signature(tf = "TranscriptionFactor"),
          function(tf) {
            return(tf@pwms)
          })

setGeneric("combine",
           function(tf1, tf2, name=name(tf1), ...) {
             if (length(list(...)) > 0) {
               combine(tf1, do.call(combine, list(tf2,...)))
             } else {
               standardGeneric("combine")
             }
           }
)
setMethod("combine",
          signature(tf1 = "TranscriptionFactor",
                    tf2 = "TranscriptionFactor"),
          function(tf1, tf2, name = NULL, ...) {
              new_name <- name
              if (is.null(name)) new_name <- name(tf1)
              pwms <- do.call('c',
                              unlist(lapply(c(tf1,tf2),
                              FUN = function(x) pwms(x))))
              return(TranscriptionFactor(name = new_name, pwms = pwms))
          }
)

setGeneric("name<-", function(tf, value) standardGeneric("name<-"))

setReplaceMethod("name",signature(tf="TranscriptionFactor",
                                  value="character"),
                 function(tf, value){
                   tf@name <- value
                   tf
                 }
)
# -----------------------------------------------------------------------------

#' @rdname Clique
#' @export
setClass("Clique",
         contains = "TranscriptionFactor",
         representation(members = "character",
                        hash = "md5")
         )

Clique <- function(tfs, name = NULL) {
  members <- sort(unlist(apply(tfs, FUN = function(x) name(x))))
  if (is.null(name)) {
    name <- paste(members,collapse = ",")
    hash <- substr(openssl::md5(name),1,8)
  }
  new_clique <- do.call(combine,tfs)
  new("Clique", new_clique, members = members, hash = hash)
}
# -----------------------------------------------------------------------------

#' @rdname CliqueList
#' @export
setClass("CliqueList",
         representation(Cliques = "character"))

#' @rdname CRCView
#' @export
setClass("CRCView",
         contains = "RangedSummarizedExperiment",
         representation(name = "character"))
