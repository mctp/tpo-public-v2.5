#' Convert dtVault to (lazy) dpVault
#'
#' @param dtv dtVault
#' @return dpVault
#'
#' @export
dpVault <- function(dtv) {
    dpv <- list()
    ## maps
    dpv$maps  <- lapply(dtv$maps, function(dtt) {
        lazy_dt(dtt)
    })
    ## tables
    dpv$tables  <- lapply(dtv$tables, function(dtt) {
        lazy_dt(dtt)
    })
    ## meta
    dpv$meta <- lazy_dt(dtv$meta)
    ## anno
    dpv$anno <- lazy_dt(dtv$anno)
    dpv$storage <- "tbl_lazy"
    dpv$format <- dtv$format
    return(dpv)
}

#' Collect dpVault to dtVault
#'
#' @param dpv dpVault
#' @param index boolean; whether to restore dtVault indexes
#' @return dtVault
#'
#' @export
dpVaultCollect <- function(dpv, index=TRUE) {
    dtv <- list()
    ## maps
    dtv$maps  <- lapply(dpv$maps, function(dtt) {
        lazy_dt(dtt) |> as.data.table()
    })
    ## tables
    dtv$tables  <- lapply(dpv$tables, function(dtt) {
        lazy_dt(dtt) |> as.data.table()
    })
    ## meta
    dtv$meta <- dpv$meta |> as.data.table()
    ## anno
    dtv$anno <- dpv$anno |> as.data.table()
    dtv$storage <- "data.table"
    dtv$format <- dpv$format
    if (index) {
        dtv <- dtVaultIndex(dtv)
    }
    return(dtv)
}
