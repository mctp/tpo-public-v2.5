.modifyDuckDBTable <- function(verb, con, dt, table_name, method="direct") {
    if (!is.null(dt)) {
        if (method=="direct") {
            # preserves factors, and supports list columns
            if (verb=="create") {
                duckdb::dbWriteTable(con, table_name, dt)
            } else if (verb=="append") {
                duckdb::dbAppendTable(con, table_name, dt)
            }
        } else if (method=="parquet") {
            ## alternative method
            if (verb=="create") {
                verb <- c("CREATE TABLE", "AS")
            } else if (verb=="append") {
                verb <- c("INSERT INTO", "")
            }
            temp_fn <- tempfile()
            write_parquet(dt, temp_fn)
            SQL <- sprintf("%s %s %s SELECT * FROM read_parquet('%s')", verb[1], table_name, verb[2], temp_fn)
            dbSendQuery(con, SQL)
            unlink(temp_fn)
        }
    }
}

.initializeDuckDBTable <- function(con, dt, table_name) {
    .modifyDuckDBTable("create", con, dt, table_name)
}

.updateDuckDBTable <- function(con, dt, table_name) {
    .modifyDuckDBTable("append", con, dt, table_name)
}

.modifyDuckDB <- function(fun, con, vault) {
    for (group in c("maps", "tables")) {
        dt_group <- vault[[group]]
        for (dt_name in names(dt_group)) {
            dt <- dt_group[[dt_name]]
            dt_name_sql <- str_replace(dt_name, "\\.", "_")
            table_name <- paste(group, dt_name_sql, sep="_")
            fun(con, dt, table_name)
        }
    }
    fun(con, vault$meta, "meta")
}

.initializeDuckDB <- function(con, vault) {
    .modifyDuckDB(.initializeDuckDBTable, con, vault)
    .initializeDuckDBTable(con, vault$anno, "anno")
}

.updateDuckDB <- function(con, vault) {
    .modifyDuckDB(.updateDuckDBTable, con, vault)
}

#' Create DuckDB from legacy vault or dtVault RDS files
#'
#' @param dbfn DuckDB filename
#' @param rds_fns Vault RDS files
#' @param anno dtVaultAnno
#' @param verbose display progress
#' @return NULL
#'
#' @export
createDuckDBFromVaultFiles <- function(dbfn, rds_fns, anno=NULL, verbose=TRUE) {
    con <- dbConnect(duckdb::duckdb(), dbdir = dbfn)
    for (i in seq_along(rds_fns)) {
        rds_fn <- rds_fns[[i]]
        if (verbose) {
            message(sprintf("loading: %s %s/%s", rds_fn, i, length(rds_fns)))
        }
        vault <- readRDS(rds_fn)
        vault <- dtVaultUpdate(vault, anno)
        if (i==1) {
            .initializeDuckDB(con, vault)
        } else {
            .updateDuckDB(con, vault)
        }
    }
    dbDisconnect(con, shutdown=TRUE)
    duckdb_shutdown(duckdb::duckdb())
    return(NULL)
}

#' Create DuckDB from a legacy vault or dtVault
#'
#' @param dbfn DuckDB filename
#' @param vault legacy or dtVault
#' @param anno dtVaultAnno object
#' @return NULL
#'
#' @export
createDuckDBFromVault <- function(dbfn, vault, anno=NULL) {
    vault <- dtVaultUpdate(vault, anno)
    con <- dbConnect(duckdb::duckdb(), dbdir = dbfn)
    .initializeDuckDB(con, vault)
    dbDisconnect(con, shutdown=TRUE)
    duckdb_shutdown(duckdb::duckdb())
    return(NULL)
}
