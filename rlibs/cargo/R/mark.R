.variantMark <- function(tbl) {
  if (!is.null(tbl)) {
    mark <- tbl[,.(
      id,
      var_id,
      var_mark=ifelse(!is.na(HGVSp), HGVSp, HGVSc)
    )]
    setkey(mark, id, var_id)
    return(mark)
  }
  return(NULL)
}

somaticMark <- function(repr) {
  .variantMark(repr)
}

ccfMark <- function(tbl) {
  if (!is.null(tbl)) {
    mark <- tbl[,.(
      id,
      var_id,
      var_mark=ifelse(!is.na(ccf), ccf, pmin(AFT*2, 1.0))
    )]
    setkey(mark, id, var_id)
    return(mark)
  }
  return(NULL)
}

germlineMark <- function(repr) {
  .variantMark(repr)
}

structuralMark <- function(repr) {
  if (!is.null(repr)) {
    mark <- rbind(
      repr[,.(id, var_id, end="e1", var_mark=CSQ1)],
      repr[,.(id, var_id, end="e2", var_mark=CSQ2)]
    )
    setkey(mark, id, var_id, end)
    return(mark)
  }
  return(NULL)
}

fusionMark <- function(repr) {
  if(!is.null(repr)) {
    mark <- rbind(
      repr[,.(id, var_id, end="f5", var_mark=CSQF5)],
      repr[,.(id, var_id, end="f3", var_mark=CSQF3)],
      repr[,.(id, var_id, end="r5", var_mark=CSQR5)],
      repr[,.(id, var_id, end="r3", var_mark=CSQR3)]
    )
    setkey(mark, id, var_id, end)
    return(mark)
  }
  return(NULL)
}

segmentMark <- function(repr) {
  if(!is.null(repr)) {
    mark <- repr[,.(
      id,
      var_id,
      var_mark=seg
    )]
    setkey(mark, id, var_id)
    return(mark)
  }
  return(NULL)
}

geneCopyMarkC <- function(tbl) {
  if(!is.null(tbl)) {
    mark <- tbl[,.(
      id,
      gene_id,
      C=Cmin
    )]
    setkey(mark, id, gene_id)
    return(mark)
  }
  return(NULL)
}

geneCopyMarkK <- function(tbl) {
  if(!is.null(tbl)) {
    mark <- tbl[,.(
      id,
      gene_id,
      K=Kmin
    )]
    setkey(mark, id, gene_id)
    return(mark)
  }
  return(NULL)
}

geneExpressionMark <- function(tbl) {
  if(!is.null(tbl)) {
    mark <- tbl[,.(
      id,
      gene_id,
      exp=round(trpkm, 2))
    ]
    setkey(mark, id, gene_id)
    return(mark)
  }
  return(NULL)
}

geneSymbolMark <- function(tbl) {
  if(!is.null(tbl)) {
    mark <- tbl[,.(
      id,
      gene_id,
      gene_name)
    ]
    setkey(mark, id, gene_id)
    return(mark)
  }
  return(NULL)
}
