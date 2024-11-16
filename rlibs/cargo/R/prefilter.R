structural.prefilter <- function(rep) {
    return(rep)
}

fusion.prefilter <- function(rep) {
    return(rep)
}

somatic.prefilter <- function(rep) {
    rep <- rep[(
        xfet<0.1 &
        AFT>AFN &
        Consequence!="intron_variant"
    )]
    return(rep)
}

somatic.to.prefilter <- function(rep) {
    rep <- rep[(
        DPT>8 &
        ADT>3 &
        Consequence!="intron_variant"
    )]
    return(rep)
}

germline.prefilter <- function(rep) {
    rep <- rep[(
        gnomad_af<0.01 & kg_af<0.01 &
        ADN>4 & AFN>0.1 &
        (IMPACT %in% c("MODERATE", "HIGH") | !is.na(clinid))
    )]
    return(rep)
}

germline.to.prefilter <- function(rep) {
    rep <- rep[(
        gnomad_af<0.01 & kg_af<0.01 &
        (IMPACT %in% c("MODERATE", "HIGH") | !is.na(clinid))
    )]
    return(rep)
}
