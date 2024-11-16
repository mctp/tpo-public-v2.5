#" @export
unifyData <- function(vault) {

  ## call maps
  m.som <- vault$maps$somatic
  m.grm <- vault$maps$germline
  m.str <- vault$maps$structural
  m.fus <- vault$maps$fusion
  m.seg <- vault$maps$segment

  ## call tables
  x.som <- somaticMark(vault$tables$somatic)
  x.ccf <- ccfMark(vault$tables$somatic)
  x.grm <- germlineMark(vault$tables$germline)
  x.str <- structuralMark(vault$tables$structural)
  x.fus <- fusionMark(vault$tables$fusion)
  x.seg <- segmentMark(vault$tables$segment)

  ## create gene master list
  v.som <- getGeneList(x.som,m.som,"som")
  v.ccf <- getGeneList(x.ccf,m.som,"ccf")
  v.grm <- getGeneList(x.grm,m.grm,"grm")
  v.str <- getGeneList(x.str,m.str,"str")
  v.fus <- getGeneList(x.fus,m.fus,"fus")
  v.seg <- getGeneListSeg(x.seg,m.seg)
  x.sym <- geneSymbolMark(vault$tables$gene.copy)
  x.abc <- geneCopyMarkC(vault$tables$gene.copy)
  x.abk <- geneCopyMarkK(vault$tables$gene.copy)
  x.exp <- geneExpressionMark(vault$tables$gene.expression)

  v.all <- Reduce(function(...) {merge(..., all=TRUE, by=c("gene_id", "id"))},
    list(v.som, v.ccf, v.grm, v.str, v.fus, v.seg, x.sym, x.abc, x.abk, x.exp)) %>%
    dplyr::relocate(id)
  return(v.all)
}

#" Get the gene lists for a map/table combination from vault object.
#"
#" @param x A table from a vault object
#" @param m A map from a vault object
#" @param name The type of table/map output from the vault object
#" @return The genes present in that map/table combination. Else NULL.
#" @examples
#"   getGeneList(x.som,m.som,"som")
#"   getGeneList(x.grm,m.grm,"grm")
getGeneList <- function(x,m,name) {
  if(!is.null(x) && !is.null(m)) {
    ret <- x[m][,.(temp=list(unique(var_mark))), c("gene_id", "id")]
    colnames(ret)[3] <- name
    setkey(ret, gene_id)
    return(ret)
  }
  return(NULL)
}

#" Get the gene lists for a map/table combination from vault object.
#"
#" @param x A segment table from a vault object
#" @param m A map from a vault object
#" @param name The type of table/map output from the vault object
#" @return The genes present in that map/table combination. Else NULL.
#" @examples
#"   getGeneListSeg(x.som,m.som)
#"   getGeneListSeg(x.grm,m.grm)
getGeneListSeg <- function(x,m)  {
  if(!is.null(x) && !is.null(m)) {
    tmp <- x[m]
    tmp <- tmp[,.(seg=list(.N)), c("gene_id", "id")]
    ret <- tmp[seg>1]
    setkey(ret, gene_id)
    return(ret)
  }
  return(NULL)
}

#" @export
unifyDataToChar <- function(vud) {
  vudc <- vud
  cols <- c("som", "grm", "str", "fus", "seg")
  foreach(col=cols) %do% {
    tmp <- sapply(vud[[col]], paste, collapse=",")
    tmp[nchar(tmp)==0] <- NA
    vudc[[col]] <- tmp
  }
  return(vudc)
}

#" @export
unifyDataFromChar <- function(vudc) {
  vud <- vudc
  cols <- c("som", "grm", "str", "fus", "seg")
  foreach(col=cols) %do% {
    tmp <- str_split(vudc[[col]], pattern=",")
    tmp[is.na(nchar(tmp))] <- list(NULL)
    vud[[col]] <- tmp
  }
  return(vud)
}
