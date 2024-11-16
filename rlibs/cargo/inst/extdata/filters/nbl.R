#### Implemented NBL filter functions
#' somatic.filter.fun
#' germline.filter.fun

#### Not implemented NBL filter functions
#' fusion.filter.fun
#' structural.filter.fun
#' segment.filter.fun


LIKELY_BENIGN <- c("Benign", "Benign/Likely_benign", "Likely_benign")
UNLIKELY_BENIGN <- c("Pathogenic/Likely_pathogenic", "Pathogenic", "Likely_pathogenic", "drug_response")




#### somatic

#' Runs somatic filters
#'
#' @param rep - table to filter
#' @return filtered table
#'
somatic.filter.fun <- function(REP, GOI=NULL){

  if( is.null(GOI) ) { GOI <- readLines(system.file("extdata/goi/goi-germline.txt", package="cargo")) }

  tn <- .somatic.filter.tn(REP)
  to <- .somatic.filter.to(REP, GOI)

  REP %>%
    dplyr::mutate(pass = ifelse(pass=='pass' & !(tn) & !is.na(ADN), 'F-custom-tn', pass)) %>%
    dplyr::mutate(pass = ifelse(pass=='pass' & !(to) & is.na(ADN), 'F-custom-to', pass)) %>%
    return(.)
}

#' somatic filter of data for Tumor Normal
#'
#' @param rep - table to filter
#' @return T/F values regarding pass
#'
.somatic.filter.tn <- function(REP) {
  fd <- REP[,(
    (tlod>5 & nlod>3 & flod>1) &
      (IMPACT %in% c("MODERATE", "HIGH") | (SYMBOL == "TERT" & (chr=="chr5" & pos==1295135 & alt=="A") | (chr=="chr5" & pos==1295113 & alt=="A"))) &
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
      mqrs>-3 &
      mqs>23 &
      ## germline
      (gnomad_af == 0 | (gnomad_af < 0.001 & AFT>0.025)) &
      (Nn < 0.010 | cosmic_cnt>100) &
      (Tn < 0.023 | cosmic_cnt>100)
  )]
  return(fd)
}

#' somatic filter of data for Tumor only
#'
#' @param rep - table to filter
#' @return T/F values regarding pass
#'
.somatic.filter.to <- function(REP, GOI) {
  fd <- REP[,(
    is.na(ADN) &
      #(SYMBOL %in% GOI) &
      (tlod>5) &
      (IMPACT %in% c("MODERATE", "HIGH") | (SYMBOL == "TERT" & (chr=="chr5" & pos==1295135 & alt=="A") | (chr=="chr5" & pos==1295113 & alt=="A"))) &
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
      (gnomad_af == 0 | (gnomad_af < 0.001 & AFT>0.025)) &
      (Nn < 0.01 | cosmic_cnt>100) &
      (Tn < 0.05 | cosmic_cnt>100)
  )]
  return(fd)
}




#### germline

#' Runs germline filters
#'
#' @param rep - table to filter
#' @return filtered table
#'
germline.filter.fun <- function(REP, GOI=NULL){

  if( is.null(GOI) ) { GOI <- readLines(system.file("extdata/goi/goi-germline.txt", package="cargo")) }

  tn <- .germline.filter.tn(REP, GOI)
  to <- .germline.filter.to(REP, GOI)

  REP %>%
    dplyr::mutate(pass = ifelse(pass=='pass' & tn==F & !is.na(ADN), 'F-custom-tn', pass)) %>%
    dplyr::mutate(pass = ifelse(pass=='pass' & to==F & is.na(ADN), 'F-custom-to', pass)) %>%
    return()

}

#' germline filter of data for Tumor Normal
#'
#' @param rep - table to filter
#' @return T/F values regarding pass
#'
.germline.filter.tn <- function(REP, GOI){

  fg <- REP[,(
    ## GOI
    (SYMBOL %in% GOI & !clinvar %in% LIKELY_BENIGN & ADT>2 & AFN>2 & AFT>0.1 & AFN>0.1) |
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
  )]
  ## deal w/ NaN
  fg = ifelse(is.nan(REP$AFT), F, fg)
  return(fg)
}

#' germline filter of data for Tumor only
#'
#' @param rep - table to filter
#' @return T/F values regarding pass
#'
.germline.filter.to <- function(REP, GOI){

  fg <- REP[,(
    (SYMBOL %in% GOI) &
      ## missing normal
      is.na(ADN) &
      ## impact
      (!clinvar %in% LIKELY_BENIGN) &
      (clinvar %in% UNLIKELY_BENIGN | gnomad_af <= 0.005 | kg_af <= 0.005) &
      (IMPACT %in% c("HIGH", "MODERATE")) &
         ## rarity
         AFT>0.2 &
         ADT>0.2
      )]
  return(fg)
}




#### Other (not implemented)

#' structural
#'
#' @param rep - table to filter
#' @return filtered table
#'
structural.filter.fun <- function(rep) {
  return(rep)
}

#' fusion
#'
#' @param rep - table to filter
#' @return filtered table
#'
fusion.filter.fun <- function(rep) {
  return(rep)
}

#' segment
#'
#' @param rep - table to filter
#' @return filtered table
#'
segment.filter.fun <- function(rep) {
  return(rep)
}

