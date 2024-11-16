#### annotation

#' annotation filter for variants
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.variant.annotation.fun <- function(rep, anno) {
    gids <- anno |> filter(anno_goi) |> pull(gene_id)
    tri_annotation <- rep |> mutate(tmp = !gene_id %in% gids) |> pull(tmp)
    return(tri_annotation)
}


#' annotation filter for somatic variants
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.somatic.annotation.fun <- triage.variant.annotation.fun


#' annotation filter for germline variants
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.germline.annotation.fun <- triage.variant.annotation.fun


#' annotation filter structural
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.structural.annotation.fun <- function(rep, anno) {
    gids <- anno |> pull(gene_id)
    tri_annotation <- rep |> mutate( tmp=!(gene_id.1 %in% gids | gene_id.2 %in% gids) ) |> pull(tmp)
    return(tri_annotation)
}


#' annotation filter fusion
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.fusion.annotation.fun <- function(rep, anno) {
    gids <- anno |> pull(gene_id)
    in_gids <- function(x) any(x %in% gids)
    if (nrow(rep)) {
      tri_annotation <- rep |>
        mutate( tmp=!(sapply(str_split(gene_ids.5.1, ":"), in_gids) |
                        sapply(str_split(gene_ids.3.1, ":"), in_gids) |
                        sapply(str_split(gene_ids.5.2, ":"), in_gids) |
                        sapply(str_split(gene_ids.3.2, ":"), in_gids))
        ) |>
        pull(tmp)
    } else {
      tri_annotation <- logical(0)
    }
    return(tri_annotation)
}




#### GOI

#' goi filter somatic
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.somatic.goi.fun <- function(rep, anno) {
  gids <- anno |> filter(somatic_goi) |> pull(gene_id)
  tri_goi <- rep |> mutate( tmp=!(gene_id %in% gids) ) |> pull(tmp)
  return(tri_goi)
}


#' goi filter germline
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.germline.goi.fun <- function(rep, anno) {
  gids <- anno |> filter(germline_goi) |> pull(gene_id)
  tri_goi <- rep |> mutate( tmp=!(gene_id %in% gids) ) |> pull(tmp)
  return(tri_goi)
}




#### artifact

#' artifact filter somatic
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.somatic.artifact.fun <- function(rep, anno) {
  gids <- anno |> filter(somatic_art) |> pull(gene_id)
  tri_art <- rep |> mutate( tmp=(gene_id %in% gids) ) |> pull(tmp)
  return(tri_art)
}




#### consequence

#' consequence filter general
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.variant.consequence.fun <- function(rep, anno) {
  tri_consequence <- rep |>
    mutate(IMPACT = ifelse(
      (
        ## TERT promoter variants
        (chr=="chr5" & pos==1295135 & alt=="A") |
          (chr=="chr5" & pos==1295113 & alt=="A") |
          ## TP53 p.Thr125= p.Glu224=
          (chr=="chr17" & pos==7675994) |
          (chr=="chr17" & pos==7674859)
        ## FOXA1 3' UTR
        ## CCND1 3' UTR
      ), "HIGH", IMPACT)) |>
    mutate(tmp = !(IMPACT %in% c("HIGH", "MODERATE"))) |> 
    pull(tmp)
  return(tri_consequence)
}


#' consequence filter somatic
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.somatic.consequence.fun <- triage.variant.consequence.fun


#' consequence filter germline
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.germline.consequence.fun <- triage.variant.consequence.fun




#### phase

#' phase filter somatic
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.somatic.phase.fun <- function(rep, anno) {
  fail <- rep |>
    mutate( key = paste0(id, '--', gene_id, '--', pid) ) |>
    filter( !is.na(pid) & duplicated(key) ) |>
    pull(key)
  tri_phase <- rep |>
    mutate( key = paste0(id, '--', gene_id, '--', pid) ) |>
    mutate( tmp = key %in% fail ) |>
    pull(tmp)
  return(tri_phase)
}




#### evidence

#' evidence filter fusion
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.fusion.evidence.fun <- function(rep, anno) {
  tri_evidence <- rep |> 
    mutate(tmp = !((hq.bpt & (mm2.valid|gmap.valid)) | hi.bpt)) |> 
    pull(tmp)
  return(tri_evidence)
}


#' evidence filter fusion
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.structural.evidence.fun <- function(rep, anno) {
  tri_evidence <- rep |> mutate(tmp = is.na(MANTA.FILTER)) |> pull(tmp)
  return(tri_evidence)
}


#' evidence filter somatic
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.somatic.evidence.fun <- function(rep, anno) {
  GOI <- anno |> filter(somatic_goi) |> pull(gene_id)
  ## Pass TN
  fdtn <- .somatic.filter.data.tn.fun(rep)
  fatn <- .somatic.filter.anno.tn.fun(rep, GOI)
  ## Pass TO
  fdto <- .somatic.filter.data.to.fun(rep)
  fato <- .somatic.filter.anno.to.fun(rep, GOI)
  ## apply
  tri_evidence <- rep |> mutate( tmp=ifelse(is.na(ADN), (fdto&fato), (fdtn&fatn)) ) |> pull(tmp) 
  return(tri_evidence)
}


#' evidence filter germline
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.germline.evidence.fun <- function(rep, anno) {
  GOI <- anno |> filter(germline_goi) |> pull(gene_id)
  ## filter TN
  ftn <- .germline.filter.tn.fun(rep)
  ## filter TO
  fto <- .germline.filter.to.fun(rep)
  ## if ADN missing pick TO
  tri_evidence <- rep |>
    mutate( tmp = ifelse(is.na(ADN), fto, ftn) ) |>
    pull(tmp)
  return(tri_evidence)
}


#' germline filter for tumor-normal
#'
#' @param rep germline dtVault table
#' @return boolean vector where TRUE == fail
#'
.germline.filter.tn.fun <- function(rep) {
  fg <- rep |>
    mutate( tmp = 
              (
                ## impact
                !(clinvar %in% LIKELY_BENIGN) &
                  (AFN >= 0.2 & ADN >= 8 & DPN >= 40) &
                  (
                    ## rare or pathogenic
                    (gnomad_af <= 0.005 & kg_af <= 0.005) | (clinvar %in% UNLIKELY_BENIGN)
                  )
              )
    ) |>
    mutate(tmp=!tmp) |>
    pull(tmp)
  return(fg)
}


#' germline filter for tumor-only
#'
#' @param rep germline dtVault table
#' @return boolean vector where TRUE == fail
#'
.germline.filter.to.fun <- function(rep) {
  fg <- rep |>
    mutate( tmp=
              (
                ## impact
                !(clinvar %in% LIKELY_BENIGN) &
                  (
                    ## rare or pathogenic
                    (gnomad_af <= 0.005 & kg_af <= 0.005) | (clinvar %in% UNLIKELY_BENIGN)
                  )
              )
    ) |>
    mutate(tmp=!tmp) |>
    pull(tmp)
  return(fg)
}


#' germline filter for genes of interest
#'
#' @param rep germline dtVault table
#' @return boolean vector where TRUE == fail
#'
.germline.filter.goi.fun <- function(rep, GOI=NULL) {
  f <- rep |> mutate( tmp=(gene_id %in% GOI) ) |> pull(tmp)
  return(f)
}


#' somatic filter for tumor-normal
#'
#' @param rep somatic dtVault table
#' @return boolean vector where TRUE == fail
#'
.somatic.filter.data.tn.fun <- function(rep) {
  fd <- rep |>
    mutate( tmp=
              ## failsafe filters, variants passing those strict filters should pass regardless
              (   ## SNVs and indels outside long STRs and RMSKs
                ((!str | str_len<4) & (!is.indel | is.na(rmsk_hit))) &
                  ## balanced failsafe filters
                  (tlod>12 & nlod>18 & (flod>9 | (ADN==0 & DPN>=100)) & xfet<1e-2 & ADT_FWD>=6 & ADT_REV>=6 & ADN<=2 & AFT>AFN*6)
              ) |
              (   ## indels inside long STRs and RMSKs
                ((str & str_len>=4) | (is.indel & !is.na(rmsk_hit))) &
                  ## strict failsafe filters
                  (tlod>24 & nlod>24 & (flod>12 | (ADN==0 & DPN>=100)) & xfet<5e-3 & ADT_FWD>=9 & ADT_REV>=9 & ADN<=1 & AFT>AFN*8)
              ) |
              ## default filters
              (   ## SNVs and indels outside long STRs and RMSKs
                ((!str | str_len<4) & (!is.indel | is.na(rmsk_hit))) &
                  ## balanced complex filters
                  ADT>4 &
                  (ADN<=3 | ((ADN / (DPN+0.1))<=0.0066)) &
                  AFT>=0.024 &
                  (AFN<0.012 | (AFN<0.02 & ADN==1) | (AFN<0.025 & ADN==0)) &
                  tlod>6 &
                  nlod>6 &
                  flod>3 &
                  xfet<0.05 &
                  (ADT_FWD>1 | (ADT_FWD==1 & ADT<10) | (ADT_FWD==1 & ADT<=12 & sor < 1.5)) &
                  (ADT_REV>1 | (ADT_REV==1 & ADT<10) | (ADT_FWD==1 & ADT<=12 & sor < 1.5)) &
                  sor<4 &
                  ((is.na(mqrs) & AFT>0.9) | (!is.na(mqrs) & mqrs>-3)) &
                  mqs>24.75
              ) |
              (
                ## indels inside long STRs and RMSKs
                ((str & str_len>=4) | (is.indel & !is.na(rmsk_hit))) &
                  ## strict complex filters
                  ADT>9 &
                  (ADN<=1 | ((ADN / (DPN+0.1))<=0.0066)) &
                  AFT>=0.045 &
                  AFN<0.005 &
                  tlod>12 &
                  nlod>12 &
                  flod>6 &
                  xfet<0.025 &
                  ADT_FWD>2 &
                  ADT_REV>2 &
                  sor<3.5 &
                  ((is.na(mqrs) & AFT>0.9) | (!is.na(mqrs) & mqrs>-3)) &
                  mqs>27
              )
    ) |>
    mutate(tmp=!tmp) |>
    pull(tmp)
  return(fd)
}


#' somatic filter for tumor-only
#'
#' @param rep somatic dtVault table
#' @return boolean vector where TRUE == fail
#'
.somatic.filter.data.to.fun <- function(rep) {
  fd <- rep |>
    mutate( tmp =
              (   ## SNVs and indels outside long STRs and RMSKs
                ((!str | str_len<4) & (!is.indel | is.na(rmsk_hit))) &
                  ## balanced filters
                  ADT>8 &
                  # ADN<=3 &
                  AFT>=0.04 &
                  (ADT_FWD>1 | (ADT_FWD==1 & ADT<10)) &
                  (ADT_REV>1 | (ADT_REV==1 & ADT<10)) &
                  sor<4 &
                  ((is.na(mqrs) & AFT>0.9) | (!is.na(mqrs) & mqrs>-3)) &
                  mqs>25
              ) |
              (
                ## indels inside long STRs and RMSKs
                ((str & str_len>=4) | (is.indel & !is.na(rmsk_hit))) &
                  ## strict filters
                  ADT>9 &
                  AFT>=0.05 &
                  ADT_FWD>2 &
                  ADT_REV>2 &
                  sor<3.5 &
                  ((is.na(mqrs) & AFT>0.9) | (!is.na(mqrs) & mqrs>-3)) &
                  mqs>27
              )
    ) |>
    mutate(tmp=!tmp) |>
    pull(tmp)
  return(fd)
}


#' somatic filter of annotation data for Tumor normal
#'
#' @param rep somatic dtVault table
#' @param GOI somatic genes of interest
#' @return boolean vector where TRUE == fail
#'
.somatic.filter.anno.tn.fun <- function(rep, GOI) {
  fa <- rep |>
    mutate( tmp=
              (
                ## important genes / variants
                (cosmic_cnt > 25 | gene_id %in% GOI | clinvar %in% UNLIKELY_BENIGN) &
                ## not a common germline variant or recurrent artifact
                (all.pon>=3 | clinvar %in% UNLIKELY_BENIGN) &
                ## near-absence in the normal
                ((ADN==0) | (ADN==1 & AFN<=0.006) | (ADN==2 & AFN<=0.004)) &
                (
                  (
                    ## low evidence outside of long STRs
                    ((!str | str_len<4) & (!is.indel | is.na(rmsk_hit))) &
                      ((ADT>2 & AFT>0.021) | (ADT>3 & AFT>0.019) | (ADT>4 & AFT>0.017) | (ADT>5 & AFT>0.015) |
                       (ADT>6 & AFT>0.013) | (ADT>7 & AFT>0.011) | (ADT>8)) &
                      (tlod>6 & nlod>6 & ADT_FWD>=1 & ADT_REV>=1) &
                      ((sblr<2.5 & sor<3.5) | (sor < 1.5)) & 
                      ((!is.na(mqrs) & mqrs>-2.5) & mqs>26.0 & ecnt<5)
                  ) |
                  (
                    ## higher evidence within long STRs
                    ((str & str_len>=4) | (is.indel & !is.na(rmsk_hit))) &
                      ((ADT>5 & AFT>0.027) | (ADT>6 & AFT>0.025) | (ADT>7 & AFT>0.023) | (ADT>8 & AFT>0.021)) &
                      (tlod>9 & nlod>9 & ADT_FWD>=2 & ADT_REV>=2) &
                      ((sblr<1.75 & sor<3.0) | (sor < 1.5)) &
                      ((!is.na(mqrs) & mqrs>-2.0) & mqs>26.0 & ecnt<3)
                  )
                )
              )
    ) |> 
    mutate(tmp=!tmp) |>
    pull(tmp)
  return(fa)
}


#' somatic filter using annotation data for tumor-only
#'
#' @param rep somatic dtVault table
#' @param GOI somatic genes of interest
#' @return boolean vector where TRUE == fail
#'
.somatic.filter.anno.to.fun <- function(rep, GOI) {
  fa <- rep |>
    mutate( tmp=
              (   ## important genes / variants
                (clinvar %in% UNLIKELY_BENIGN | cosmic_cnt > 25 | gene_id %in% GOI) &
                  ## SNVs and indels outside long STRs and RMSKs
                  ((!str | str_len<5) & (!is.indel | is.na(rmsk_hit))) &
                  ## low-evidence but strict quality and absence in normal
                  ADT>2 &
                  AFT>0.02 &
                  ADT_FWD>=1 &
                  ADT_REV>=1 &
                  sor<2.5 &
                  (!is.na(mqrs) & mqrs>-0.5) &
                  mqs>26 &
                  (gnomad_af==0 & kg_af == 0 & is.na(pon))
              )
    ) |>
    mutate(tmp=!tmp) |>
    pull(tmp)
  return(fa)
}




#### custom

#' Not implemented but needed to run
#'
#' @param rep somatic dtVault table
#' @param anno dxVault annotation
not.implemented <- function(rep, anno){
  to_return <- F
  return(to_return)
}


triage.somatic.custom.fun <- not.implemented
triage.germline.custom.fun <- not.implemented
triage.structural.custom.fun <- not.implemented
triage.fusion.custom.fun <- not.implemented



