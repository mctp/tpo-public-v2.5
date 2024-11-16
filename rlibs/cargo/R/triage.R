#' Determines filter status based on annotation, evidence, consequence, artifact, and pid
#'
#' @param dxv dbVault dpVault dtVault
#' @param anno custom dtVaultAnno object; optional
#' @param fset custom filter set (filename or name); optional
#' @param update boolean; update vault tables [default: TRUE]
#' @param somatic boolean; whether to triage somatic
#' @param germline boolean; whether to triage germline
#' @param structural boolean; whether to triage structural
#' @param fusion boolean; whether to triage fusion
#'
#' @return dxVault or triage result
#'
#' @export
vaultTriage <- function(dxv, anno=NULL, fset=NULL, somatic=TRUE, germline=TRUE,
  structural=TRUE, fusion=TRUE, update=TRUE) {
  ## if provided anno overrides vault anno
  if (is.null(anno)) {
    anno <- dxv$anno
  }
  ## import default filter functions
  source(system.file("extdata/triage/default.R", package="cargo"))
  ## if provided import filtering functions which will override defaults
  if (!is.null(fset)) {
    if (file.exists(fset)) {
      source(fset)
    } else {
      source(system.file(sprintf("extdata/triage/%s.R", fset), package="cargo"))
    }
  }
  ## run triaging
  tri_somatic <- tri_germline <- tri_structural <- tri_fusion <- NULL
  if (somatic) {
    tri_somatic <- .triageSomatic(dxv$tables$somatic, dxv$anno)
  }
  if (germline) {
    tri_germline <- .triageGermline(dxv$tables$germline, dxv$anno)
  }
  if (structural) {
    tri_structural <- .triageStructural(dxv$tables$structural, dxv$anno)
  }
  if (fusion) {
    tri_fusion <- .triageFusion(dxv$tables$fusion, dxv$anno)
  }
  if (update) {
    ## update vault tables
    dxv$tables$somatic <- .updateTableFilter(dxv$tables$somatic, tri_somatic)
    dxv$tables$germline <- .updateTableFilter(dxv$tables$germline, tri_germline)
    dxv$tables$structural <- .updateTableFilter(dxv$tables$structural, tri_structural)
    dxv$tables$fusion <- .updateTableFilter(dxv$tables$fusion, tri_fusion)
    return(dxv)
  } else {
    ## return triaging result
    triage <- list(somatic=tri_somatic, germline=tri_germline,
                   structural=tri_structural, fusion=tri_fusion)
    return(triage)
  }
}

#' Removes variants based on triaging result.
#'
#' @param dxv dbVault dpVault dtVault
#' @param rules variant filtering rules
#' @return dxVault
#'
#' default rules:
#'   list(
#'     somatic="annotation|consequence|evidence|artifact",
#'     germline="annotation|consequence|evidence|goi",
#'     structural="evidence",
#'     fusion="evidence"
#'   )
#'
#' @export
vaultFilter <- function(dxv, rules=NULL) {
  if (is.null(rules)) {
    rules <- list(
      somatic="annotation|consequence|evidence|artifact|manual|custom",
      germline="annotation|consequence|evidence|goi|manual|custom",
      structural="evidence|manual|custom",
      fusion="evidence|manual|custom"
    )
  }
  if (!is.list(rules)){
    rules <- split(unname(rules),names(rules))
  }
  for (table_name in names(rules)) {
    rule <- rules[[table_name]]
    tab <- dxv$tables[[table_name]]
    map <- dxv$maps[[table_name]]
    if (dxv$storage == "data.table") {
      tabf <- tab[ !grepl(rule, triage) | grepl('keep', triage) ]
      mapf <- map[tabf[,.(id, var_id)], on=c("id", "var_id")]
    } else {
      tabf <- tab |> filter( !grepl(rule, triage) | grepl('keep', triage) )
      mapf <- map |> semi_join(tabf, by=c("id", "var_id"))
    }
    dxv$tables[[table_name]] <- tabf
    dxv$maps[[table_name]] <- mapf
  }
  return(dxv)
}

#' Updates table filter column based on triaging result.
#'
#' Note: for DB-backed tables update is always in-place and
#' expects non-subsetted dbVault
#'
#' @param rep variant table
#' @param tri triaging table
#' @param in_place update data.table in-place?
#'
#' @return rep
.updateTableFilter <- function(rep, tri, in_place=TRUE) {
  
  ## get existing triage vals
  existing <- rep$triage
  existing[is.na(existing)] <- ''
  
  if (!is.null(tri)) {
    ## create filter string
    id_vars <- c("id", "var_id")
    tri_vars <- setdiff(colnames(tri), id_vars)
    tri_name <- list()
    for (tri_var in tri_vars) {
        tri_name[[tri_var]] <- ifelse(tri[[tri_var]], tri_var, NA_character_)
    }
    pf <- tri_name |> data.frame() |> unite(".", sep=";", na.rm=TRUE) |> pull() %>% paste0(existing, ';', .) %>% gsub('^;|;$|;;', '', .)
    #pf[pf==""] <- "PASS"
    tri$triage <- pf
    ## update filtering string
    tri_upd <- tri[,c("id", "var_id", "triage")]
    if ("tbl_duckdb_connection" %in% class(rep)) {
      stopifnot(in_place)
      ## update database in place
      con <- rep[[1]][[1]]
      duckdb::duckdb_register(con, "tri", tri_upd)
      rep |> rows_update(tbl(con, "tri"), by=c("id", "var_id"), unmatched="ignore", in_place=TRUE)
      duckdb::duckdb_unregister(con, "tri")
    } else if ("data.table" %in% class(rep) && in_place) {
      ## update data.table in place
      tri_upd <- as.data.table(tri_upd)
      rep[tri_upd, triage:=i.triage, on=c("id", "var_id")]
    } else if ("dtplyr_step" %in% class(rep) && in_place) {
      ## update data.table in place when using dpVault
      rep <- rep |> select(-triage) |> left_join(tri_upd, by=c("id", "var_id")) |> relocate(triage)
    } else {
      ## not in-place
      rep <- rep |> rows_update(tri_upd, by=c("id", "var_id"))
    }
  }
  return(rep)
}

#' Determines triage status of somatic variants
#'
#' @param somatic somatic table
#' @param anno anno table
#'
#' @return triaged table
.triageSomatic <- function(somatic, anno) {
  ## convert dbVault table to dtVault if necessary
  if ("tbl_dbi" %in% class(somatic)) {
    warning(paste("collecting", somatic[[2]]$x[[1]]))
    somatic <- somatic |> collect() |> data.table()
    anno <- anno |> collect() |> data.table()
  }
  ## triage functions
  annotation <- triage.somatic.annotation.fun(somatic, anno)
  evidence <- triage.somatic.evidence.fun(somatic, anno)
  consequence <- triage.somatic.consequence.fun(somatic, anno)
  artifact <- triage.somatic.artifact.fun(somatic, anno)
  goi <- triage.somatic.goi.fun(somatic, anno)
  phase <- triage.somatic.phase.fun(somatic, anno)
  custom <- triage.somatic.custom.fun(somatic, anno)
  ## triage result table
  tri_somatic <- somatic |> select(id, var_id) |> collect() |>
    mutate(annotation, evidence, consequence, goi, artifact, phase, custom)
  if (anyNA(tri_somatic)) {
    print(which(!complete.cases(tri_somatic)))
    stop()
  }
  return(tri_somatic)
}

#' Determines triage status of germline variants
#'
#' @param germline germline table
#' @param anno anno table
#'
#' @return filtered vault
.triageGermline <- function(germline, anno) {
  if ("tbl_dbi" %in% class(germline)) {
    warning(paste("collecting", germline[[2]]$x[[1]]))
    germline <- germline |> collect() |> data.table()
    anno <- anno |> collect() |> data.table()
  }
  ## triage functions
  annotation <- triage.germline.annotation.fun(germline, anno)
  evidence <- triage.germline.evidence.fun(germline, anno)
  consequence <- triage.germline.consequence.fun(germline, anno)
  goi <- triage.germline.goi.fun(germline, anno)
  custom <- triage.germline.custom.fun(germline, anno)
  ## triage result table
  tri_germline <- germline |> select(id, var_id) |> collect() |>
    mutate(annotation, evidence, consequence, goi)
  stopifnot(!anyNA(tri_germline))
  return(tri_germline)
}

#' Determines triage status of structural variants
#'
#' @param structural structural table
#' @param anno anno table
#'
#' @return triaged table
.triageStructural <- function(structural, anno) {
  if ("tbl_dbi" %in% class(structural)) {
    warning(paste("collecting", structural[[2]]$x[[1]]))
    structural <- structural |> collect() |> data.table()
    anno <- anno |> collect() |> data.table()
  }
  ## triage functions
  annotation <- triage.structural.annotation.fun(structural, anno)
  evidence <- triage.structural.evidence.fun(structural, anno)
  custom <- triage.structural.custom.fun(structural, anno)
  ## triage result table
  tri_structural <- structural |> select(id, var_id) |> collect() |>
    mutate(annotation, evidence, custom)
  stopifnot(!anyNA(tri_structural))
  return(tri_structural)
}

#' Determines triage status of gene fusions
#'
#' @param fusion fusion table
#' @param anno anno table
#'
#' @return triaged table
.triageFusion <- function(fusion, anno) {
  if ("tbl_dbi" %in% class(fusion)) {
    warning(paste("collecting", fusion[[2]]$x[[1]]))
    fusion <- fusion |> collect() |> data.table()
    anno <- anno |> collect() |> data.table()
  }
  ## triage functions
  annotation <- triage.fusion.annotation.fun(fusion, anno)
  evidence <- triage.fusion.evidence.fun(fusion, anno)
  custom <- triage.fusion.custom.fun(fusion, anno)
  ## triage result table
  ## bug in duckdb mutate, so collect before
  tri_fusion  <- fusion |> select(id, var_id) |> collect() |>
    mutate(annotation, evidence)
  stopifnot(!anyNA(tri_fusion))
  return(tri_fusion)
}




#### Rescue variant functions

#' Changes triage value of variants based on provided table
#'
#' @param V - vault
#' @param TBL a table containing id,var_id,table,& pass value
#' @return updated vault
#' @export
dxRescueVariants <- function(V, TBL) {
  
  ## Error handling
  if( is.null(TBL) ){ stop("No variant table provided; consider using .changeVariantPass for single variant") }
  if( !all(c('triage','id','var_id','table') %in% colnames(TBL)) ){ stop("The following cols needed in TBL: triage, id, var_id, table") }
  if( !all(TBL$table %in% names(V$tables)) ){ stop("A provided table is not present in vault") }
  notInVault <- .dxRescueCheck(V, TBL)
  if( length(notInVault)>0 ){ stop( paste0("The following id--var_id pairs are not in the vault:  ", paste0( notInVault, collapse = ', ')) )}
  ## select
  som_vals <- TBL %>% filter(table=='somatic')
  grm_vals <- TBL %>% filter(table=='germline')
  str_vals <- TBL %>% filter(table=='structural')
  fus_vals <- TBL %>% filter(table=='fusion')
  ## rescue
  som_table = V$tables$somatic %>% .dxRescueVariants(som_vals, .)
  grm_table = V$tables$germline %>% .dxRescueVariants(grm_vals, .)
  str_table = V$tables$structural %>% .dxRescueVariants(str_vals, .)
  fus_table = V$tables$fusion %>% .dxRescueVariants(fus_vals, .)
  ## assign
  V$tables$somatic <- som_table
  V$tables$germline <- grm_table
  V$tables$structural <- str_table
  V$tables$fusion <- fus_table
  
  return(V)
}

#' Changes pass/fail value of a vault variant table based on provided values
#'
#' @param NEW_VALS - the variant table with new values
#' @param VAULT_TABLE the vault table needing to be updated w/ new values
#' @return updated vault
#' @export
.dxRescueVariants <- function(NEW_VALS, VAULT_TABLE){
  nms <- colnames(VAULT_TABLE)
  tbl <- NEW_VALS %>%
    ## add/reformat
    merge(., VAULT_TABLE, by=c('id', 'var_id'), all.y=T, ) %>%
    mutate(triage = paste(triage.x, triage.y, sep=';')) %>%
    mutate(triage = gsub('^NA;', '', triage)) %>%
    select(all_of(nms))
  return(tbl)
}

#' Returns id--var_id pairs from teh table that are not present in the vault
#'
#' @param V - vault
#' @param TBL a table containing id,var_id,table,& pass value
#' @return vector of id-var_id pairs not present in the vault but in the table
.dxRescueCheck <- function(V, TBL){
  ## get keys
  som = V$tables$somatic |> mutate(key = paste0(id,'--',var_id)) |> pull(key)
  grm = V$tables$germline |> mutate(key = paste0(id,'--',var_id)) |> pull(key)
  str = V$tables$structural |> mutate(key = paste0(id,'--',var_id)) |> pull(key)
  fus = V$tables$fusion |> mutate(key = paste0(id,'--',var_id)) |> pull(key)
  all = c(som,grm,str,fus)
  keys <- TBL |> mutate(key = paste0(id,'--',var_id)) |> pull(key)
  ## check if in vault
  missing = keys[!keys %in% all]
  return(missing)
}
