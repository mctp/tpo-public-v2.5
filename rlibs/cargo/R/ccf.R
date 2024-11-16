#' Calculates (meta)vault CCF.
#'
#' @param DTV dtVault
#' @param CFG cargo config via cargoConfig()
#' @param FIELD which table to calculate for; default somatic
#' @return vault with CCF/M in selected table
#'
#' @export
calculateVaultCCF <- function(DTV, CFG, FIELD="somatic") {
  dxVaultCheck(DTV, storage = "data.table", format="v2")
  ids <- DTV$meta |> pull(id)
  mv <- lapply(ids, dxVaultSubsetById, dxv=DTV) %>%
    lapply(., .calculateVaultCCF, CFG, FIELD) %>%
    dtVaultStack(.)
  return(mv)
  }

#### Some relevant math
# (eq.1) ccf = vaf / (purity * m) * (purity * C + (1 - purity) * nC)
# (eq.2) vaf = ccf * purity * m / (purity * C + (1 - purity) * nC)
# (eq.3) purity = 1 / (((m * ccf / vaf) - C) / nC + 1)
# Quick way of computing CCF:
# u_i <- ((ADT / DPT) / purity) * (purity * C + (1 - purity) * nC)
# m_i <- pmax(1, round(u_i))
# ccf_i <- u_i / m_i

#' Calculates vault CCF for an individual vault
#'
#' @param V The vault
#' @param CFG cargo config via cargoConfig()
#' @param FIELD which table to calculate for; default somatic
#' @return vault with CCF/M in selected table
#'
.calculateVaultCCF <- function(V, CFG, FIELD="somatic") {
  tbl_field <- getVaultTable(V, FIELD)
  tbl_seg <- getVaultTable(V, "segment")
  ## handle null variant table/no variants
  stopifnot(!is.null(tbl_field))
  stopifnot(!is.null(tbl_seg))
  if (nrow(tbl_field) == 0 || nrow(tbl_seg) == 0) {
    tbl_field[,ccf:=numeric()]
    tbl_field[,m:=numeric()]
    v <- setVaultTable(V, FIELD, tbl_field)
    return(v)
  }
  ## prepare data
  ccfd <- .prepareVaultCCF(V, FIELD)
  dl <- list(snv = ccfd, purity = V$meta$purity, ploidy = V$meta$ploidy)
  m <- .calcM(dl, CFG)
  ccf <- (dl$snv$AFT * ((dl$snv$nC * (1 - dl$purity)) + (dl$snv$C * dl$purity))) / (dl$purity * m)
  ccf <- ifelse(is.infinite(ccf), NA_real_, ccf)
  ## add to vault
  tmp <- getVaultTable(V, FIELD)
  tmp$ccf <- ccf
  tmp$m <- m
  setkey(tmp, id, var_id)
  v <- setVaultTable(V, FIELD, tmp)
  ## return
  return(v)
}

#' Prepares mutation table for CCF calculation. Robust to MV
#'
#' @param V The vault
#' @param FIELD which table
#' @return prepared data table
#'
.prepareVaultCCF <- function(V, FIELD) {
  ## mutations
  tbl <- getVaultTable(V, FIELD)
  ## purity/ploidy
  tmp <- distinct(V$meta[,c("id", "purity", "ploidy")])
  tbl <- merge(tbl, tmp, by="id", all.x=TRUE)
  ## pos
  tbl$start <- tbl$pos
  tbl$end   <- tbl$pos
  ## segment
  segment <- getVaultTable(V, "segment")
  setkey(segment, id, chr, start, end)
  ## overlap
  hits <- foverlaps(tbl, segment, by.x = c("id", "chr", "start", "end"), by.y = c("id", "chr", "start", "end"), nomatch = NULL, which = TRUE)
  tbl[["C"]] <- NA_real_
  tbl[["K"]] <- NA_real_
  tbl[["lr"]] <- NA_real_
  tbl[["nlr"]] <- NA_real_
  tbl[["sC"]] <- NA_real_
  tbl[["len"]] <- NA_real_
  tbl[["nC"]] <- NA_real_
  tbl[["seg"]] <- NA_real_
  tbl[hits$xid,"C"] <- segment[hits$yid,"C"]
  tbl[hits$xid,"K"] <- segment[hits$yid,"K"]
  tbl[hits$xid,"lr"] <- segment[hits$yid,"lr"]
  tbl[hits$xid,"nlr"] <- segment[hits$yid,"nlr"]
  tbl[hits$xid,"sC"] <- segment[hits$yid,"sC"]
  tbl[hits$xid,"len"] <- segment[hits$yid,"len"]
  tbl[hits$xid,"nC"] <- segment[hits$yid,"nC"]
  tbl[hits$xid,"seg"] <- segment[hits$yid,"seg"]
  ## return
  return(tbl)
}

#' Calculates multiplicity of variant
#'
#' @param DL data list
#' @param CFG cargo config via cargoConfig()
#' @return vector of Ms - multiplicities
#'
.calcM <- function(DL, CFG) {
  ## compute expected VAFs for all possible M's
  mc <- .ccfMGrid(DL)
  af <- as.data.table(DL$snv[, c("chr", "start", "end", "AFT", "DPT", "seg", "C", "K")])
  af$idx <- seq_len(nrow(af))
  colnames(af)[colnames(af) == "C"] <- "correctC"
  colnames(af)[colnames(af) == "K"] <- "correctK"
  ##
  tmp1 <- rbindlist(rep(list(af), nrow(mc)))[order(idx)]
  tmp2 <- rbindlist(rep(list(mc), nrow(af)))
  afc <- cbind(tmp1, tmp2)
  ##
  afc[, `:=`(Ef, Ef.aut)]
  afc[chr %in% c("chrX", "chrY"), `:=`(Ef, Ef.sex)]
  afc[, `:=`(c("Ef.aut", "Ef.sex"), NULL)]
  afc <- afc[C == correctC, ]
  afc[, `:=`(likelihood = dbeta(Ef,
                                shape1 = AFT * pmin(DPT, CFG$ccf.dp.af.max) + 1,
                                shape2 = (1 - AFT) * pmin(DPT, CFG$ccf.dp.af.max) + 1,
                                log = TRUE), likelihood.type = "Beta")]
  afc <- afc[order(plyr::desc(likelihood)), .SD[1], keyby = .(idx)]
  af[afc, on = "idx", `:=`(M, i.M)]
  M <- af$M
  ## return
  return(M)
}

#' function required to calculate multiplicity
#'
#' @param DL data list
#'
.ccfMGrid <- function(DL) {
  purity <- DL$purity
  C <- DL$snv$C
  sex <- ifelse(1 %in% unique(DL$snv$nC), "male", "female")
  ## make a grid
  max.C <- max(C,na.rm = TRUE)
  MC <- as.data.table(expand.grid(M=1:max.C,C=1:max.C))[M<=C]
  if (sex == "female") {
    MC[,":="(
      Ef.aut=(purity * M) / (purity * C + 2 * (1-purity)),
      Ef.sex=(purity * M) / (purity * C + 2 * (1-purity))
    )]
  } else {
    MC[,":="(
      Ef.aut=(purity * M) / (purity * C + 2 * (1-purity)),
      Ef.sex=(purity * M) / (purity * C + 1 * (1-purity)) # male has 1 copy on normal sex chr
    )]
  }
  return(MC)
}


.allelicModelData <- function(v, table, varcol, segcol) {
  segr <- with(v$tables$segment, GRanges(
    chr, IRanges(start, end), "*"
  ))
  hascol <- colnames(v$tables$segment)
  if (!all(segcol %in% hascol)) {
    warning(sprintf("Some segment columns are missing %s", paste(setdiff(segcol, hascol), collapse=", ")))
    segcol <- intersect(segcol, hascol)
  }
  mcols(segr) <- v$tables$segment[,segcol, with=FALSE]
  varr <- with(v$tables[[table]], GRanges(
    chr, IRanges(pos, pos), "*"
  ))
  hascol <- colnames(v$tables[[table]])
  if (!all(varcol %in% hascol)) {
    warning(sprintf("Some variant columns are missing %s", paste(setdiff(varcol, hascol), collapse=", ")))
    varcol <- intersect(varcol, hascol)
  }
  mcols(varr) <- v$tables[[table]][,varcol, with=FALSE]
  ## find nearest segment for each somatic variant
  mcols(varr) <- cbind(mcols(varr), mcols(segr[nearest(varr, segr)]))
  mcols(varr)$purity <- v$meta$purity
  mcols(varr)$ploidy <- v$meta$ploidy
  allcol <- c("seqnames", varcol, segcol, "purity", "ploidy")
  var <- as.data.table(varr)[,allcol, with=FALSE]
  setnames(var, "seqnames", "chr")
  return(var)
}

#' Estimate purity from variants and CNV model
#'
#' @param dtv dtVault (1 sample)
#' @param strategy named heuristic used to estimate purity [diploid-modal]
#' @return purity (numeric)
#'
#' @export
estimatePurityFromSomaticVariants <- function(dtv, strategy="diploid-modal") {
  dxVaultCheck(dtv, storage="data.table")
  dxVaultCheckCardinality(dtv, 1)
  table <- "somatic"
  varcol <- c("AFT")
  segcol <- c("C", "K", "sC", "nC")
  adm <- .allelicModelData(dtv, table, varcol, segcol)
  ## filter sub-clonal
  if (strategy=="diploid-modal") {
    ## pre-filter to diploid regions
    admf <- adm[nC==2]
    aftf <- adm$AFT
    ## test if data is bimodal between 0 and 0.5 VAF
    hasdip <- dip.test(aftf[aftf < 0.5])
    if (hasdip$p.value < 0.05) {
      ## remove subclonal variants
      aftf <- aftf[aftf > hasdip$statistic]
    }
    modaft <- mlv(aftf, method="shorth")
    m <- 1
    ccf <- 1
    C <- 2
    nC <- 2
    purity <- 1/((((m*ccf/modaft)-C)/nC)+1)
  } else {
    stop("Strategy unknown.")
  }
  return(purity)
}







#' Estimate purity from variants and CNV model using MVF as input
#'
#' @param MVF metavault (1 sample)
#' @param strategy named heuristic used to estimate purity [diploid-modal]
#' @return table (ID, numeric purity)
#'
#' @export
estimatePurityfromMetavault <- function(MVF, strategy="diploid-modal", priority_genes = NULL) {
  
  # MVF$tables$somatic <- MVF$tables$somatic %>% filter(pass == "pass") %>% data.frame()
  
  tmp.somatic <- MVF$tables$somatic %>% filter(pass == "pass") %>% data.frame()
  tmp.segment <- MVF$tables$segment %>% data.frame()
  tmp.meta <- MVF$meta
  table <- "somatic"
  varcol <- c("AFT", "m", "ccf", "id")
  segcol <- c("C", "K", "sC", "nC", "id")
  
  
  mvfallelicModelData <- function(tmp.somatic,tmp.segment,tmp.meta, varcol, segcol, pID = NULL) {
    
    if(!is.null(pID)){tmp.segment <- tmp.segment %>% filter(id == pID) %>% data.frame()}
    
    segr <- with(tmp.segment, GRanges(
      chr, IRanges(start, end), "*"
    ))
    hascol <- colnames(tmp.segment)
    if (!all(segcol %in% hascol)) {
      warning(sprintf("Some segment columns are missing %s", paste(setdiff(segcol, hascol), collapse=", ")))
      segcol <- intersect(segcol, hascol)
    }
    mcols(segr) <- tmp.segment[,segcol]
    varr.df <- data.frame(tmp.somatic)
    if(!is.null(pID)){varr.df <- varr.df %>% filter(id == pID)%>% data.frame()}
    varr <- with(varr.df, GRanges(
      chr, IRanges(pos, pos), "*"
    ))
    hascol <- colnames(tmp.somatic)
    if (!all(varcol %in% hascol)) {
      warning(sprintf("Some variant columns are missing %s", paste(setdiff(varcol, hascol), collapse=", ")))
      varcol <- intersect(varcol, hascol)
    }
    varr.mcol <- data.frame(tmp.somatic[,c("id", varcol)])
    if(!is.null(pID)){varr.mcol <- varr.mcol %>% filter(id == pID)%>% data.frame()}
    mcols(varr) <- varr.mcol
    ## find nearest segment for each somatic variant
    mcols(varr) <- cbind(mcols(varr), mcols(segr[nearest(varr, segr)]))
    mcols(varr)$purity <- ifelse(is.null(pID), tmp.meta$purity, tmp.meta$purity[tmp.meta$id == pID])
    mcols(varr)$ploidy <- ifelse(is.null(pID), tmp.meta$ploidy, tmp.meta$ploidy[tmp.meta$id == pID])
    allcol <- unique(c("seqnames", varcol, segcol, "purity", "ploidy"))
    var <- as.data.table(varr)
    setnames(var, "seqnames", "chr")
    return(var)
  }
  adm <- mvfallelicModelData(tmp.somatic,tmp.segment,tmp.meta, varcol, segcol)
  get_purity <- function(adm, strategy){
    ## filter sub-clonal
    
    if (strategy=="diploid-modal") {
      ## pre-filter to diploid regions
      admf <- adm %>% filter(nC==2) %>% data.frame()
      aftf <- adm$AFT
      ## test if data is bimodal between 0 and 0.5 VAF
      hasdip <- dip.test(aftf[aftf < 0.6])
      if (hasdip$p.value < 0.05) {
        ## remove subclonal variants
        aftf <- aftf[aftf > hasdip$statistic]
      } else {
        aftf <- aftf[aftf >= mean(aftf)]
      } 
      if(length(aftf) == 0){return(NA)}
      modaft <- mlv(aftf, method="shorth")
      m <- 1
      ccf <- 1
      C <- 2
      nC <- 2
      purity <- 1/((((m*ccf/modaft)-C)/nC)+1)
    } else {
      stop("Strategy unknown.")
    }
    return(purity)
  }
  
  if(length(unique(MVF$meta$id)) > 1){
    print("more than one sample detected")
  }
  purity_table <- data.frame()
  for(i in unique(MVF$meta$id)){
    p <- get_purity(data.frame(adm %>% filter(id == i)), strategy)
    purity_table <- rbind(purity_table, data.frame(ID = i, purity = p))
  }
  
  return(purity_table)
  
}