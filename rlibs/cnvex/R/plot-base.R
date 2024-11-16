CNVEX_COPY_COL <- c(
    K0="#78b09c", KN="#000000", CN="#690033",
    C0="#4500AC", C1="#6B58EE", C2="#000000", C3="#FC9272", C4="#FB6A4A",
    C5="#EF3B2C", C6="#CB181D", C7="#A50F15", C8="#67000D"
)

CNVEX_SEG_COL <- c("#6495ED", "#DD8080", "#CDE2B8")

segCol <- function(seg) {
    out <- list(palette=CNVEX_SEG_COL, keys=factor((seg %% 3)+1))
    return(out)
}

## C,K,lr to color
copyCol <- function(C=NULL, K=NULL, lr=NULL, wt.col=NULL, copy.col=NULL) {
    if (is.null(copy.col)) {
        copy.col <- CNVEX_COPY_COL
    }
    if (!is.null(wt.col)) {
        copy.col["C2"] <- wt.col
        copy.col["KN"] <- wt.col
    }
    if (!is.null(C) && !is.null(K)) {
        K[C==1] <- 0
        key <- ifelse(K %in% 0 & C>=2, "K0", ifelse(C<=8, paste0("C", C), "CN"))
    }
    if (!is.null(C) && is.null(K)) {
        key <- ifelse(C<=8, paste0("C", C), "C8")
    }
    if (is.null(C) && !is.null(K)) {
        key <- ifelse(K %in% 0, "K0", "KN")
    }
    if (!is.null(lr)) {
        cr <- colorRampPalette(c(copy.col["C0"], copy.col["C2"], copy.col["C8"]))
        copy.col <- cr(99)
        lr[lr < -2] <- -2
        lr[lr > +2] <- +2
        key <- cut(lr, breaks = seq(-2, +2, len=100), include.lowest=TRUE)
    }
    out <- list(palette=copy.col, keys=key)
    return(out)
}

colTypeSwitch <- function(obj, col.type, wt.col, decode=TRUE) {
    if (col.type=="seg") {
        out <- segCol(obj$seg)
    } else if (col.type=="C") {
        out <- copyCol(C=as.integer(obj$C), wt.col=wt.col)
    } else if (col.type %in% c("CK", "KC")) {
        out <- copyCol(C=as.integer(obj$C), K=as.integer(obj$K), wt.col=wt.col)
    } else if (col.type=="K") {
        out <- copyCol(K=as.integer(obj$K), wt.col=wt.col)
    } else if (col.type=="lr") {
        out <- copyCol(lr=as.numeric(obj$lr), wt.col=wt.col)
    } else if (col.type=="none") {
        out <- copyCol(C=as.integer(obj$C), wt.col=wt.col)
        out$keys <- rep("C2",length(out$keys))
    } else {
        stop("Unknown col.type")
    }
    if (decode) {
        out <- out$palette[out$keys]
    }
    return(out)
}

copyToRatio <- function(C, purity, ploidy) {
    p <- purity
    P <- ploidy
    D <- (P * p) + 2 * (1 - p)
    D <- (P * p) + 2 * (1 - p)
    lr <- log2((p * C + (1 - p) * 2) / D)
    return(lr)
}

modelGrid <- function(purity, ploidy, abs.range=24) {
    p <- purity
    D <- ploidy * p + 2 * (1 - p)
    C <- seq(0, abs.range, by = 1)
    grid <- as.data.table(expand.grid(M = 0:abs.range, C = 0:abs.range))[M <= C]
    grid[, lr := log2((p * C + (1 - p) * 2) / D)]
    grid[, af := (p * M + 1 * (1 - p)) / (p * C + 2 * (1 - p))]
    return(grid[])
}

modelLines <- function(purity, Msel=0, Csel=1) {
    p <- purity
    baf <- (p * Msel + 1 * (1-p)) / (p * Csel + 2 * (1-p))
    baf.lines <- data.table(
        baf=c(baf, 1-baf),
        lab.mc=rep(paste(Msel, Csel, sep="/"), 2),
        lab.cdiff=sprintf("%+d", Csel-2)
    )
    return(baf.lines)
}
