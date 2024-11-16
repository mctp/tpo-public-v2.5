
LIKELY_BENIGN <- c("Benign", "Benign/Likely_benign", "Likely_benign")
UNLIKELY_BENIGN <- c("Pathogenic/Likely_pathogenic", "Pathogenic", "Likely_pathogenic", "drug_response")




#### not implemented

#' Not implemented but needed to run
#'
#' @param rep somatic dtVault table
#' @param anno dxVault annotation
not.implemented <- function(rep, anno){
  to_return <- F
  return(to_return)
}

# somatic
triage.somatic.evidence.fun <- not.implemented
triage.somatic.annotation.fun <- not.implemented
triage.somatic.consequence.fun <- not.implemented
triage.somatic.artifact.fun <- not.implemented
triage.somatic.goi.fun <- not.implemented
triage.somatic.phase.fun <- not.implemented

# germline
triage.germline.evidence.fun <- not.implemented
triage.germline.annotation.fun <- not.implemented
triage.germline.consequence.fun <- not.implemented
triage.germline.goi.fun <- not.implemented

# fusion
triage.fusion.evidence.fun <- not.implemented

# structural
triage.structural.evidence.fun <- not.implemented




#### somatic

#' evidence filter somatic
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.somatic.custom.fun <- function(rep, anno) {
  GOI <- anno |> filter(somatic_goi) |> pull(gene_id)
  ## Pass TN
  fdtn <- .triage.somatic.custom.tn.fun(rep)
  ## Pass TO
  fato <- .triage.somatic.custom.to.fun(rep, GOI)
  ## apply
  tri_evidence <- rep |> mutate( tmp=ifelse(is.na(ADN), (fato), (fdtn)) ) |> pull(tmp)
  return(tri_evidence)
}


#' somatic filter for tumor-normal
#'
#' @param rep somatic dtVault table
#' @return boolean vector where TRUE == fail
#'
.triage.somatic.custom.tn.fun <- function(rep) {
  fd <- rep |>
    mutate( tmp=
              (
                (tlod>5 & nlod>3 & flod>1) &
                  (IMPACT %in% c("MODERATE", "HIGH") | (gene_name == "TERT" & (chr=="chr5" & pos==1295135 & alt=="A") | (chr=="chr5" & pos==1295113 & alt=="A"))) &
                  ## hard threshold
                  (AFT>=0.02 & AFN<0.02) &
                  (ccf>=0.05 | is.na(ccf)) &
                  (ADT>=7 & ADN<7) &
                  ## evidence
                  tlod > 6 &
                  nlod > 6 &
                  flod > 1.5 &
                  ## increase
                  (xfet < 0.05 | (xfet < 0.1 & ADN==0 & tlod>7.5 & nlod>7.5)) &
                  (AFT > AFN * 10) &
                  ## artifacts
                  ADT_FWD>0 &
                  ADT_REV>0 &
                  sor<4 &
                  (!is.na(mqrs) & mqrs>-3.0) & 
                  mqs>23 &
                  ## germline
                  (gnomad_af == 0 | (gnomad_af < 0.001 & AFT>0.025))
                  #cosmic_cnt>100
                  #(Nn < 0.010 | cosmic_cnt>100) &
                  #(Tn < 0.023 | cosmic_cnt>100)
              )
    ) |>
    mutate(tmp=!tmp) |>
    pull(tmp)
  return(fd)
}


#' somatic filter for tumor-only
#'
#' @param rep somatic dtVault table
#' @param GOI somatic genes of interest
#' @return boolean vector where TRUE == fail
#'
.triage.somatic.custom.to.fun <- function(rep, GOI) {
  fa <- rep |>
    mutate( tmp=
              (
                (gene_id %in% GOI) &
                is.na(ADN) &
                  (tlod>5) &
                  (IMPACT %in% c("MODERATE", "HIGH") | (gene_name == "TERT" & (chr=="chr5" & pos==1295135 & alt=="A") | (chr=="chr5" & pos==1295113 & alt=="A"))) &
                  ## hard threshold
                  (AFT>=0.02) &
                  (ccf>=0.05 | is.na(ccf)) &
                  (ADT>=10) &
                  ## artifacts
                  ADT_FWD>0 &
                  ADT_REV>0 &
                  sor<4 &
                  ( (!is.na(mqrs)) & mqrs>-3) &
                  ( (!is.na(mqs)) & mqs>23) &
                  ## germline
                  (gnomad_af == 0 | (gnomad_af < 0.001 & AFT>0.025))
                  #cosmic_cnt>100
                  #(Nn < 0.01 | cosmic_cnt>100) &
                  #(Tn < 0.05 | cosmic_cnt>100)
              )
    ) |>
    mutate(tmp=!tmp) |>
    pull(tmp)
  return(fa)
}




#### germline

#' evidence filter germline
#'
#' @param rep dxVault table
#' @param anno dxVault annotation
#' @return boolean vector where TRUE == fail
#'
triage.germline.custom.fun <- function(rep, anno) {
  GOI <- anno |> filter(germline_goi) |> pull(gene_id)
  ## filter TN
  ftn <- .triage.germline.custom.tn.fun(rep, GOI)
  ## filter TO
  fto <- .triage.germline.custom.to.fun(rep, GOI)
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
.triage.germline.custom.tn.fun <- function(rep, GOI) {
  fg <- rep |>
  mutate( tmp = 
      ## GOI
      (gene_name %in% GOI & !clinvar %in% LIKELY_BENIGN & ADT>2 & AFN>2 & AFT>0.1 & AFN>0.1) |
        ## non-GOI
        (ADT>8 & AFT>0.2) &
        (AFN >= 0.2 & ADN >= 8 & DPN >= 40) &
        ## impact
        (!clinvar %in% LIKELY_BENIGN) &
        (clinvar %in% UNLIKELY_BENIGN | gnomad_af <= 4e-5 | kg_af <= 4e-5) &
        (IMPACT %in% c("HIGH", "MODERATE") &
           ## artifacts
           sor<4 &
           ( (gnomad_af < 0.005 & kg_af <= 0.005) | (clinvar %in% UNLIKELY_BENIGN) )
        )
    ) |>
    ## deal w/ NaN
    mutate(tmp = ifelse(is.nan(AFT), F, tmp)) |>
    mutate(tmp=!tmp) |>
    pull(tmp)
  return(fg)
}


#' germline filter for tumor-only
#'
#' @param rep germline dtVault table
#' @return boolean vector where TRUE == fail
#'
.triage.germline.custom.to.fun <- function(rep, GOI) {
  fg <- rep |>
    mutate( tmp=
              (gene_name %in% GOI) &
              ## missing normal
              is.na(ADN) &
              ## impact
              (!clinvar %in% LIKELY_BENIGN) &
              (clinvar %in% UNLIKELY_BENIGN | gnomad_af <= 0.005 | kg_af <= 0.005) &
              (IMPACT %in% c("HIGH", "MODERATE")) &
              ## rarity
              AFT>0.2 &
              ADT>0.2
    ) |>
    mutate(tmp=!tmp) |>
    pull(tmp)
  return(fg)
}



