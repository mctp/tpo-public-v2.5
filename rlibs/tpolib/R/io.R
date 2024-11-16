#' @export
list.files.null <- function(..., n=Inf) {
    ret <- list.files(...)
    if (length(ret)==0) {
        ret <- NULL
    }
    ret <- head(ret, n)
    return(ret)
}

#' @export
get.opts <- function(settings, opts) {
    ENV <- new.env(parent = .BaseNamespaceEnv)
    source(settings, local=ENV)
    opts <- merge_list(opts, ENV$OPTS)
    return(opts)
}
