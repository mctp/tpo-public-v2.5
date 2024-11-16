#' Create dbVault object from a DB connection
#'
#' @param con Supported DB connection
#' @return dbVault
#'
#' @export
dbVault <- function(con) {
    tbl_names <- dbListTables(con)
    dbv <- list(
        maps=list(),
        tables=list(),
        meta=tbl(con, "meta"),
        anno=tbl(con, "anno"),
        storage=class(con),
        format="v2"
    )
    for (tbl_name in grep("maps|tables", tbl_names, value=TRUE)) {
        tn <- str_split(tbl_name, "_", 2)[[1]]
        tn2 <- str_replace(tn[2], "_", ".")
        dbv[[tn[1]]][[tn2]] <- tbl(con, tbl_name)
    }
    return(dbv)
}

#' Collect on-disk dbVault to in-memory dtVault
#'
#' @param dbv dbVault
#' @param index boolean; whether to restore dtVault indexes
#' @return dtVault
#'
#' @export
dbVaultCollect <- function(dbv, index=TRUE, lazy=TRUE) {
    dtv <- list()
    ## maps
    dtv$maps  <- lapply(dbv$maps, function(dbt) {
        dbt |> collect() |> as.data.table()
    })
    ## tables
    dtv$tables  <- lapply(dbv$tables, function(dbt) {
        dbt |> collect() |> as.data.table()
    })
    ## meta
    dtv$meta <- dbv$meta |> collect() |> as.data.table()
    ## anno
    dtv$anno <- dbv$anno |> collect() |> as.data.table()
    dtv$storage <- "data.table"
    dtv$format <- dbv$format
    if (index) {
        dtv <- dtVaultIndex(dtv)
    }
    return(dtv)
}