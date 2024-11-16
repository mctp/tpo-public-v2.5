#' Default filter function set:
#' 
#' somatic.filter.fun
#' germline.filter.fun
#' germline.filter.tn.fun
#' germline.filter.to.fun
#' consequence.filter.fun
#' artifact.filter.fun
#' fusion.filter.fun
#' structural.filter.fun
#' segment.filter.fun
#' phase.filter.fun




#' somatic filter of data for Tumor Normal from gatk
#' 
#' @param rep - table to filter
#' @return T/F values regarding pass
#' 
.somatic.filter.data.tn.fun <- function(rep) {
  fd <- rep[,
            ## failsafe filters, variants passing those strict filters should pass regardless
            (   ## SNVs and indels outside long STRs and RMSKs
              ((!str | str_len<4) & (!is.indel | lengths(rmsk_hit)==0)) &
                ## balanced failsafe filters
                (tlod>12 & nlod>18 & xfet<1e-2 & ADN<=2 & AFT>AFN*6)
            ) |
              (   ## indels inside long STRs and RMSKs
                ((str & str_len>=4) | (is.indel & lengths(rmsk_hit)!=0)) &
                  ## strict failsafe filters
                  (tlod>24 & nlod>24 & xfet<5e-3 & ADN<=1 & AFT>AFN*8)
              ) |
              ## default filters
              (   ## SNVs and indels outside long STRs and RMSKs
                ((!str | str_len<4) & (!is.indel | lengths(rmsk_hit)==0)) &
                  ## balanced complex filters
                  ADT>4 &
                  (ADN<=3 | ((ADN / (DPN+0.1))<=0.0066)) &
                  AFT>=0.024 &
                  (AFN<0.012 | (AFN<0.02 & ADN==1)) &
                  tlod>6 &
                  nlod>6 &
                  xfet<0.05 &
                  (is.na(mqrs) & AFT>0.9)
              ) |
              (
                ## indels inside long STRs and RMSKs
                ((str & str_len>=4) | (is.indel & lengths(rmsk_hit)!=0)) &
                  ## strict complex filters
                  ADT>9 &
                  (ADN<=1 | ((ADN/ (DPN+0.1))<=0.0066)) &
                  AFT>=0.045 &
                  AFN<0.005 &
                  tlod>12 &
                  nlod>12 &
                  xfet<0.025 &
                  (is.na(mqrs) & AFT>0.9)
              )
  ]
  return(fd)
}

#' somatic filter of annotation data for Tumor normal from gatk
#' 
#' @param rep - table to filter
#' @return T/F values regarding pass
#' 
.somatic.filter.anno.tn.fun <- function(rep) {
  tier1 <- readLines(system.file("extdata/goi/goi-somatic-tier1.txt", package="cargo"))
  fa <- rep[,
            (
              ## important genes / variants
              (cosmic_cnt > 25 | SYMBOL %in% tier1 | clinvar %in% UNLIKELY_BENIGN) &
                ## not a common germline variant or recurrent artifact
                (all.pon>=3 | clinvar %in% UNLIKELY_BENIGN) &
                ## near-absence in the normal
                ((ADN==0) | (ADN==1 & AFN<=0.006) | (ADN==2 & AFN<=0.004)) &
                (
                  (
                    ## low evidence outside of long STRs
                    ((!str | str_len<4) & (!is.indel | lengths(rmsk_hit)==0)) &
                      ((ADT>2 & AFT>0.021) | (ADT>3 & AFT>0.019) | (ADT>4 & AFT>0.017) | (ADT>5 & AFT>0.015) |
                         (ADT>6 & AFT>0.013) | (ADT>7 & AFT>0.011) | (ADT>8)) &
                      (tlod>6 & nlod>6) &
                      ecnt<5
                  ) |
                    (
                      ## higher evidence within long STRs
                      ((str & str_len>=4) | (is.indel & lengths(rmsk_hit)!=0)) &
                        ((ADT>5 & AFT>0.027) | (ADT>6 & AFT>0.025) | (ADT>7 & AFT>0.023) | (ADT>8 & AFT>0.021)) &
                        (tlod>9 & nlod>9) &
                        ecnt<3
                    )
                )
            )
  ]
  return(fa)
}

#' Runs somatic filters
#' 
#' @param rep - table to filter
#' @return filtered table
#' 
somatic.filter.fun <- function(rep) {
  ## modify
  rep %<>%
    dplyr::mutate(is.indel = str_length(ref)!=str_length(alt)) %>%
    dplyr::mutate(all.pon = ((gnomad_af==0) + (lengths(pon)==0) + (kg_af==0)))
  ## Pass TN
  fdtn <- .somatic.filter.data.tn.fun(rep)
  fatn <- .somatic.filter.anno.tn.fun(rep)
  ## Pass TO
  fdto <- .somatic.filter.data.to.fun(rep)
  fato <- .somatic.filter.anno.to.fun(rep)
  ## apply
  somatic <- rep %>%
    dplyr::mutate(pass = ifelse(pass=='pass' & !(fdtn|fatn) & !is.na(ADN), 'F-evidence_tn', pass)) %>%
    dplyr::mutate(pass = ifelse(pass=='pass' & !(fdto|fato) & is.na(ADN), 'F-evidence_to', pass)) %>%
    dplyr::select(-c(is.indel,all.pon))
  return(somatic)
}





#### Germline

#' Runs germline filters
#' 
#' @param rep - table to filter
#' @return filtered table
#' 
germline.filter.fun <- function(rep) {
  ## filter TN
  ftn <- germline.filter.tn.fun(rep)
  ## filter TO
  fto <- germline.filter.to.fun(rep)
  ## filter goi
  fgoi <- germline.filter.goi(rep)
  ## return
  germline <- rep %>%
    dplyr::mutate(pass = ifelse(pass=='pass' & ftn==FALSE & !is.na(ADN), 'F-evidence_data', pass)) %>%
    dplyr::mutate(pass = ifelse(pass=='pass' & fto==FALSE &  is.na(ADN), 'F-evidence_data', pass)) %>%
    ## goi
    dplyr::mutate(pass = ifelse(pass=='pass' & fgoi==FALSE, 'F-goi', pass))
  return(germline)
}

#' germline filter for tumor-normal
#' 
#' @param rep - table to filter
#' @return T/F values regarding pass
#' 
germline.filter.tn.fun <- function(rep) {
  fg <- rep[,(
    ## impact
    !(clinvar %in% LIKELY_BENIGN) &
      (AFN >= 0.2 & ADN >= 8 & DPN >= 40) &
      (
        ## rare or pathogenic
        (gnomad_af <= 0.005 & kg_af <= 0.005) | (clinvar %in% UNLIKELY_BENIGN)
      )
  )]
  return(fg)
}

#' germline filter for Tumor only
#' 
#' @param rep - table to filter
#' @return T/F values regarding pass
#' 
germline.filter.to.fun <- function(rep) {
  fg <- rep[,(
    ## impact
    !(clinvar %in% LIKELY_BENIGN) &
      (
        ## rare or pathogenic
        (gnomad_af <= 0.005 & kg_af <= 0.005) | (clinvar %in% UNLIKELY_BENIGN)
      )
  )]
  return(fg)
}

#' germline filter for genes of interest
#' 
#' @param rep - table to filter
#' @return T/F values regarding pass
#' 
germline.filter.goi <- function(rep){
  goi <- readLines(system.file("extdata/goi/goi-germline.txt", package="cargo"))
  f <- rep[,SYMBOL %in% goi]
  return(f)
}




#### Consequence/Artifact

#' consequence filter 
#' 
#' @param rep - table to filter
#' @return filtered table
#' 
consequence.filter.fun <- function(rep) {
  rep %>%
    dplyr::mutate(IMPACT = ifelse(
      (
        ## TERT promoter variants
        (chr=="chr5" & pos==1295135 & alt=="A") |  
          (chr=="chr5" & pos==1295113 & alt=="A") |
          ## TP53 p.Thr125= p.Glu224=
          (chr=="chr17" & pos==7675994) |
          (chr=="chr17" & pos==7674859)
        ## FOXA1 3' UTR
        ## CCND1 3' UTR
      ), 'HIGH', IMPACT)) %>%
    dplyr::mutate(pass = ifelse(pass=='pass' & !IMPACT %in% c("HIGH", "MODERATE"), 'F-consequence', pass)) %>%
    return()
}

#' artifact filter 
#' 
#' @param rep - table to filter
#' @return filtered table
#' 
artifact.filter.fun <- function(rep) {
  art <- readLines(system.file("extdata/goi/art-somatic.txt", package="cargo"))
  fa <- rep %>%
    dplyr::mutate(pass = ifelse(pass=='pass' & (SYMBOL %in% art | ecnt>=10), 'F-artifact', pass))
  return(fa)
}




#### other

#' fusion filter
#' 
#' @param rep - table to filter
#' @return filtered table
#' 
fusion.filter.fun <- function(rep, goi) {
  rep %>%
    dplyr::mutate(pass = ifelse( pass=='pass' & 
                                   !( (hq.bpt & (mm2.valid|gmap.valid) ) | hi.bpt ), 
                                 'F-evidence', pass)) %>%
    return()
}

#' structural filter 
#' 
#' @param rep - table to filter
#' @return filtered table
#' 
structural.filter.fun <- function(rep) {
  rep %>% 
    dplyr::mutate(pass = ifelse(pass=='pass' & is.na(MANTA.FILTER), 'F-evidence', pass)) %>%
    return()
}

#' segment filter 
#' 
#' @param rep - table to filter
#' @return filtered table
#' 
segment.filter.fun <- function(rep) {
  ## TODO: is there a meaningful way to filter segments here?
  return(rep)
}

#' Removes variants on same allele of gene after ordering by quality
#' Quality currently based on IMPACT/tlod
#' 
#' @param rep - table to filter
#' @return filtered table
#' 
phase.filter.fun <- function(rep) {
  ## temp variable
  rep$id_varid = paste0(rep$id, '--', rep$var_id)
  
  ## identify duplicates
  to_fail = rep %>%
    ## tmp vars
    dplyr::mutate(tmp_impact = factor(IMPACT, levels = c('LOW','MODIFIER','MODERATE','HIGH'))) %>%
    dplyr::mutate(id_pid = paste0(id, '--', pid)) %>%
    ## filter
    dplyr::filter(pass=='pass' & !is.na(pid)) %>%
    dplyr::arrange(desc(tmp_impact), desc(tlod)) %>%
    dplyr::filter(duplicated(id_pid))
  
  ## mark/return
  rep %>%
    dplyr::mutate(pass = ifelse(id_varid %in% to_fail$id_varid, 'F-pid', pass)) %>%
    dplyr::select(-id_varid) %>%
    return()
  
}
