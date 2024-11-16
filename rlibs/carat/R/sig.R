.sig.pass96 <- function(rep) {
    pass <- rep[,
                !is.na(sbs96) &
                ## artifacts
                ADT_FWD>=2 &
                ADT_REV>=2 &
                sor<4 &
                mqrs>-3 &
                mqs>23 &
                ## not in the germline
                ((ADN == 0 & nlod>6) | nlod>20) &
                Nn == 0 &
                lengths(rep$pon)==0 &
                ## present in the tumor
                AFT > 0.015 &
                tlod > 6
                ]
    return(pass)
}

sigCount <- function(rep, sig96) {
    cts96 <- rep[,.N,sbs96]
    setkey(cts96, sbs96)
    cts96 <- cts96[J(rownames(sig96))]
    cts96[is.na(N),N:=0]
    sig <- cts96[,N]
    names(sig) <- cts96$sbs96
    return(sig)
}
