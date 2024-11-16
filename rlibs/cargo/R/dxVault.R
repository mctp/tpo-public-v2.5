#' Subset dxVault by column and value
#'
#' If a column is not present in a table it is a no-op
#'
#' @param dxv dbVault dpVault dtVault
#' @param filter.col column to filter on
#' @param filter.val values to keep
#' @return dxVault
#'
.dxVaultSubset <- function(dxv, filter.col, filter.val) {
    out <- dxv
    ## maps
    out$maps  <- lapply(dxv$maps, function(dxt) {
        if (filter.col %in% names(dxt)) {
            dxt |> filter(!!as.symbol(filter.col) %in% filter.val)
        } else {
            dxt
        }
    })
    ## tables
    out$tables  <- lapply(dxv$tables, function(dxt) {
        if (filter.col %in% names(dxt)) {
            dxt |> filter(!!as.symbol(filter.col) %in% filter.val)
        } else {
            dxt
        }
    })
    ## meta
    if (filter.col %in% names(dxv$meta)) {
        out$meta <- dxv$meta |> filter(!!as.symbol(filter.col) %in% filter.val)
    }
    return(out)
}

#' Subset dxVault by id
#'
#' @param dxv dbVault dpVault dtVault
#' @param ids one or multipe ids
#' @return dxVault
#'
#' @export
dxVaultSubsetById <- function(dxv, ids) {
    filter.col <- "id"
    filter.val <- ids
    .dxVaultSubset(dxv, filter.col, filter.val)
}

#' Subset dxVault by index
#'
#' @param dxv dbVault dpVault dtVault
#' @param idxs one or multipe indexes
#' @return dxVault
#'
#' @export
dxVaultSubsetByIndex <- function(dxv, idxs) {
    filter.col <- "id"
    filter.val <- dtv$meta[idxs, id]
    .dxVaultSubset(dxv, filter.col, filter.val)
}


#' Subset dxVault by gene
#'
#' @param dxv dbVault dpVault dtVault
#' @param gene_names one or multiple gene symbols
#' @param gene_ids one or multipe gene ids
#' @return dxVault
#'
#' @export
dxVaultSubsetByGene <- function(dxv, gene_names=NULL, gene_ids=NULL) {
    if (!is.null(gene_names)) {
        gene_ids2 <- dxv$anno |> filter(gene_name %in% gene_names) |> collect() |> pull(gene_id)
    }
    genes <- c(gene_ids, gene_ids2)
    out <- dxv
    for (tbl_name in names(dxv$map)[!sapply(dxv$map, is.null)]) {
        map_sub <- dxv$maps[[tbl_name]] |> filter(gene_id %in% genes)
        tbl_sub <- dxv$tables[[tbl_name]] |> semi_join(map_sub, by=c("id", "var_id"))
        out$maps[[tbl_name]] <- map_sub
        out$tables[[tbl_name]] <- tbl_sub
    }
    for (tbl_name in c("gene.expression", "gene.copy")) {
        out$tables[[tbl_name]] <- dxv$tables[[tbl_name]] |> filter(gene_id %in% genes)
    }
    return(out)
}

#' Subset dxVault by table
#'
#' @param dxv dbVault dpVault dtVault
#' @param tables one or multipe table names
#' @return dxVault
#'
#' @export
dxVaultSubsetByTable <- function(dxv, tables) {
    out <- dxv
    to_remove <- setdiff(names(dxv$tables), tables)
    for (tbl_name in to_remove) {
        out$tables[[tbl_name]] <- out$tables[[tbl_name]] |> filter(FALSE)
        if (!(tbl_name %in% c("gene.expression", "gene.copy"))) {
            out$maps[[tbl_name]] <- out$maps[[tbl_name]] |> filter(FALSE)
        }
    }
    return(out)
}

#' Apply function to dxVault subsets
#'
#' @param dxv dbVault dpVault dtVault
#' @param fun function operating on vault
#' @param split subset by 'id' or 'gene'
#' @return dxVault
#'
#' @export
dxVaultApply <- function(dxv, fun, split="id") {
    lapply(dxv$meta$ids, function(id) {
    })
}

#' Collect on-disk dbVault or in-memory dpVault
#'
#' @param dbv dbVault
#' @param ... pass-through parameters to dbVaultCollect and dpVaultCollect
#' @return dtVault
#'
#' @export
dxVaultCollect <- function(dxv, ...) {
    if (dxv$storage %in% c("tbl_lazy")) {
        dtv <- dpVaultCollect(dxv, ...)
    } else if (dxv$storage %in% c("duckdb")) {
        dtv <- dbVaultCollect(dxv, ...)
    } else {
        dtv <- dxv
    }
    return(dtv)
}

#' Check if vault has allowed storage and version
#'
#' By default allows all backends
#'
#' @param dxv dbVault dpVault dtVault
#' @param storage allowed storage type [data.table,duckdb,tbl_lazy]
#' @param format allowed format version [v2]
#' @return NULL
#'
#' @export
dxVaultCheck <- function(dxv, storage=c("data.table", "duckdb", "tbl_lazy"), format=c("v2")) {
    if (!(dxv$storage %in% storage)) {
        stop(sprintf("Unsupported storage: %s, valid: %s", dxv$storage, paste(storage, collapse=", ")))
    }
    if (!(dxv$format %in% format)) {
        stop(sprintf("Unsupported storage: %s, valid: %s", dxv$format, paste(format, collapse=", ")))
    }
}

#' Check if vault has allowed storage and version
#'
#' @param dxv dbVault dpVault dtVault
#' @param expect expected cardinality default: 1
#' @return NULL
#'
#' @export
dxVaultCheckCardinality <- function(dxv, expect=1) {
    nsamp <- dxv$meta |> tally() |> pull(n)
    if (nsamp != expect) {
        stop(sprintf("Incorrect cardinality: %s, expect: %s", nsamp, expect))
    }
}
