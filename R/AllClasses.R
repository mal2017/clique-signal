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
           function(tf1, tf2, name=NULL, ...) {
             if (length(list(...)) > 0) {
               combine(x, do.call(combine, list(tf1,...)))
             } else {
               standardGeneric("combine")
             }
           }
)

setMethod("combine",
          signature(tf1 = "TranscriptionFactor",
                    tf2 = "TranscriptionFactor"),
          function(tf1,tf2, name = NULL, ...) {
            if (name(tf1) != name(tf2)) {
              if (is.null(name)) stop("Supply a name for your combined TF.")
              new_name <- name
            } else {
              if (is.null(name)) {
                new_name <- name(tf1)
              } else {
                new_name <- name
              }
            }
            pwms <- do.call('c',
                            unlist(lapply(c(tf1,tf2),
                                          FUN = function(x) pwms(x))))
            return(TranscriptionFactor(name = name, pwms = pwms))
})
# -----------------------------------------------------------------------------

#' @rdname Clique
#' @export
setClass("Clique",
         contains = "SimpleList",
         representation(members = "character",
                        name = "character"))


Clique <- function(tfs, name = NULL) {
  members <- sort(unlist(apply(tfs, FUN = function(x) name(x))))
  name <-
  new("Clique", SimpleList(tfs),members )

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
