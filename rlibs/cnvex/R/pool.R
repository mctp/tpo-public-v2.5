.poolFilter <- function(cov, nvar, target, hq, hi.zero, pool.filter, hi.nvar) {
    f <- logical(length(target))
    if (pool.filter == FALSE) {
        f[target] <- TRUE
    } else {
        # zero filter
        f[target] <- rowMeans(cov[target,] == 0, na.rm = TRUE) < hi.zero & rowMads(cov[target,], na.rm = TRUE)!=0
        # nvar filter
        f[target] <- f[target] & (nvar[target] < hi.nvar)
        # hq filter
        f[target] <- f[target] & (hq[target])
    }

    # handle NAs
    f[is.na(f)] <- FALSE

    return(f)
}

.outlierMaskSel <- function(xs, sd) {
    ## identify outlier tiles
    xs.med <- rowMedians(xs, na.rm = TRUE)
    xs.mad <- rowMads(xs, na.rm = TRUE)
    xs.out <- (xs - xs.med) / xs.mad
    ## lo-threshold
    xs[xs.out < -sd] <- NA_real_
    xs.min <- rowMins(xs, na.rm = TRUE)
    xs <- apply(xs, 2, function(col) {
        col[is.na(col)] <- xs.min[is.na(col)]
        return(col)
    })
    ## hi-threshold
    xs[xs.out > sd] <- NA_real_
    xs.max <- rowMaxs(xs, na.rm = TRUE)
    xs <- apply(xs, 2, function(col) {
        col[is.na(col)] <- xs.max[is.na(col)]
        return(col)
    })
    return(xs)
}

.outlierSexChrAdjust <- function(cov, sex, pool.sex, sex.chr) {
    cov[sex.chr,pool.sex != sex] <- NA_real_
    return(cov)
}

.outlierMask <- function(cov, filter, sd, sex, pool.sex, sex.chr) {
    cov[filter,] <- .outlierMaskSel(cov[filter,], sd)
    cov[!filter,] <- NA_real_
    cov <- .outlierSexChrAdjust(cov,sex,pool.sex,sex.chr)
    return(cov)
}

.fixLogRatioBounds <- function(xsr) {
    xsrm <- max(min(xsr[is.finite(xsr)]), -4)
    xsrx <- min(max(xsr[is.finite(xsr)]),  4)
    xsr[is.nan(xsr)] <- 0    ## 0 / 0
    xsr[xsr < xsrm] <-  xsrm ## 0 / x
    xsr[xsr > xsrx] <-  xsrx ## x / 0
    xsr[is.na(xsr)] <- 0
    return(xsr)
}

.medianLogRatioSel <- function(xs) {
    xsm <- rowMedians(xs,na.rm = TRUE)
    xsr <- log2(xs/xsm)
    ## fix bounds
    xsr <- .fixLogRatioBounds(xsr)
    return(xsr)
}

.medianLogRatio <- function(x, s, f, hq, gc, chr.list, opts) {
    x[(s & f),] <- .medianLogRatioSel(x[(s & f),])
    x[(!s & f),] <- .medianLogRatioSel(x[(!s & f),])
    x.biasnorm <- foreach(i = 1:ncol(x)) %dopar% {
        return(.correctBias(x[,i],s,hq,"gc",gc,opts))
    }
    x <- do.call(cbind,x.biasnorm)
    # x[f,] <- t(t(x[f,]) - colMedians(x[f,], na.rm = TRUE))
    x[f,] <- .colNormalization(x[f,],chr.list[f],opts)
    return(x)
}

.colNormalization <- function(x, chr.list,opts) {
    colnorms <- lapply(X = unique(chr.list), function(X) {
        return(t(t(x[chr.list == X,]) - colMedians(as.matrix(x[chr.list == X,]), na.rm = TRUE)))
    })
    x.norm <- do.call(rbind,colnorms)
    return(as.matrix(x.norm))
}

.medianCoverage <- function(cov) {
    m <- rowMedians(cov, na.rm = TRUE)
    return(m)
}

.importPoolData <- function(cnv.fns, gobj, opts) {
    imported_pool <- foreach(file_num = 1:length(cnv.fns)) %dopar% {
        cnv <- readRDS(cnv.fns[file_num])
        cov <- cnv$tile$n.cov
        sex <- cnv$sex
        sex.chr <- as.character(seqnames(cnv$tile)) %in% gobj$sexchr
        target <- cnv$tile$target
        hq <- cnv$tile$hq
        chr <- as.character(seqnames(cnv$tile))
        gc <- cnv$tile$gc
        reptime <- cnv$tile$reptime
        return(list(cov, sex, sex.chr, target, hq, chr, gc, reptime))
    }
    cov <- list()
    sex <- list()
    sex.chr <- list()
    target <- list()
    hq <- list()
    chr <- list()
    gc <- list()
    reptime <- list()
    for (i in 1:length(imported_pool)){
        l <- imported_pool[[i]]
        cov[[i]] <- l[[1]]
        sex[[i]] <- l[[2]]
    }

    cov <- matrix(unlist(cov), ncol = length(cov), byrow = FALSE)
    sex.chr <- l[[3]]   #these measures are the same among all samples
    target <- l[[4]]    #these measures are the same among all samples
    hq <- l[[5]]        #these measures are the same among all samples
    chr <- l[[6]]       #these measures are the same among all samples
    gc <- l[[7]]        #these measures are the same among all samples
    reptime <- l[[8]]   #these measures are the same among all samples
    pd <- list(cov=cov, sex=unlist(sex), sex.chr=sex.chr, target=target, hq=hq, chr.list=chr, gc=gc, reptime=reptime)
    pd <- .separatePoolDataBySex(pd, opts)
    return(pd)
}

.separatePoolDataBySex <- function(pd, opts) {
    ## seperate information
    cov.autosomal <- pd$cov[pd$sex.chr == FALSE,]
    cov.male <- pd$cov[pd$sex.chr == TRUE,]
    cov.male[, pd$sex == "female"] <- NA_real_
    cov.male <- rbind(cov.autosomal,cov.male)
    cov.female <- pd$cov[pd$sex.chr == TRUE,]
    cov.female[, pd$sex == "male"] <- NA_real_
    cov.female <- rbind(cov.autosomal, cov.female)

    ## join new info
    pd.new <- list()
    pd.new$male$cov <- cov.male
    pd.new$female$cov <- cov.female
    pd.new$sex.chr <- pd$sex.chr
    pd.new$target <- pd$target
    pd.new$hq <- pd$hq
    pd.new$sex <- pd$sex
    pd.new$chr.list <- pd$chr.list
    pd.new$gc <- pd$gc
    pd.new$reptime <- pd$reptime

    return(pd.new)
}

.calcPoolProjection <- function(lr.aut,lr.sex,opts) {
    if (opts$pool.method %in% c("pca","pca-nosex")) {
        p.aut <- svd(lr.aut)$u
        p.sex <- svd(lr.sex)$u
    } else if (opts$pool.method %in% c("ica", "ica-nosex")) {
        p.aut <- fastICA(lr.aut, opts$pool.n.comp, method="C")$S
        p.sex <- fastICA(lr.sex, opts$pool.n.comp, method="C")$S
    } else {
        p.aut <- NULL
        p.sex <- NULL
    }
    return(list(p.aut = p.aut, p.sex = p.sex))
}

.removeCoverageBias <- function(pd, opts) {
    cov.male <- pd$male$cov[,pd$sex == "male"]
    cov.female <- pd$female$cov[,pd$sex == "female"]

    ## in case the cohort has only one sex
    adjust.coef <- ifelse(length(unique(pd$sex))>1, .calcCovBiasCoef(cov.male,cov.female,pd$sex.chr,pd$hq,pd$target,opts), 1)

    pd$female$cov <- pd$female$cov*adjust.coef
    pd$adjust.coef <- adjust.coef
    return(pd)
}

.calcCovBiasCoef <- function(cov.male, cov.female, sex.chr, hq, target, opts) {
    cov.male.auto <- mean(cov.male[!sex.chr & hq & target,])
    cov.female.auto <- mean(cov.female[!sex.chr & hq & target,])
    return(cov.male.auto/cov.female.auto)
}

.createPool <- function(pd, opts) {
    ## 1. Remove coverage bias on autosomes between XX and XY
    ##   - compute ratio of male/female autosome coverage
    ##   - multiply femal coverage to match male
    ##   - no coverage bias if pd has only one sex
    pd <- .removeCovefrageBias(pd, opts)
    ## 2. Calculate per-tile variance (sd) actually
    ##   TODO: refactor to use rowSds from matrixStats
    nvar <- .calcNvar(pd, opts)
    ## 3. Compute threshold of hi-variance, pool.hi.nvar is a percentile e.g. 0.95
    ##   TODO: this will fail in female-only pools
    hi.nvar <- quantile(nvar$male[!pd$sex.chr], opts$pool.hi.nvar) ## thr only on autosomes, same for both sex

    ## female
    ## 4. Creates boolean filter to only use high-quality in statistical pool
    ##  - only active if pool.filter is TRUE
    ##  - even if pool.filter is FALSE removes off-target tiles
    ##  - pick tiles which are rarely 0, low-variance, and hq
    ##  TODO: figure out why filter can contain NAs
    filter.female <- .poolFilter(pd$female$cov, nvar$female, pd$target, pd$hq, opts$pool.hi.zero, opts$pool.filter, hi.nvar)
    ## 5. Makes multiple adjustments to coverage for it to work well in statistical pool
    ##  - within filter-pass tile coverages remove outliers based on pool.sd.out
    ##  - TODO: replace filter-fail tile coverages with NAs
    ##  - in a female pool mask chrX and chrY coverages in males (and vice-versa)
    cov.female <- .outlierMask(pd$female$cov, filter.female, opts$pool.sd.out, "female", pd$sex, pd$sex.chr)
    ## 6. Calculates median coverage AFTER adjustments, this means NAs for all filter-fail coverages.
    mc.female <- .medianCoverage(cov.female)
    lr.female <- .medianLogRatio(cov.female, pd$target, filter.female, pd$hq, pd$gc, pd$chr.list, opts)

    ## male
    filter.male <- .poolFilter(pd$male$cov, nvar$male, pd$target, pd$hq, opts$pool.hi.zero, opts$pool.filter, hi.nvar)
    cov.male <- .outlierMask(pd$male$cov, filter.male, opts$pool.sd.out, "male", pd$sex, pd$sex.chr)
    mc.male <- .medianCoverage(cov.male)
    lr.male <- .medianLogRatio(cov.male, pd$target, filter.male, pd$hq, pd$gc, pd$chr.list, opts)

    ## seperate autosomes & autosomes+sex.chr
    lr.female.aut <- lr.female[filter.female & !pd$sex.chr,]
    lr.female.sex <- lr.female[filter.female, pd$sex == "female"]

    lr.male.aut <- lr.male[filter.male & !pd$sex.chr,]
    lr.male.sex <- lr.male[filter.male, pd$sex == "male"]

    ## projection
    projection <- list(p.female.aut=NULL,p.female.sex=NULL,p.male.aut=NULL,p.male.sex=NULL)
    if (all(pd$sex == "female")) {
        ssprojection <- .calcPoolProjection(lr.female.aut,lr.female.sex,opts)
        projection$p.female.aut <- ssprojection$p.aut
        projection$p.female.sex <- ssprojection$p.sex
    } else if (all(pd$sex == "male")) {
        ssprojection <- .calcPoolProjection(lr.male.aut,lr.male.sex,opts)
        projection$p.male.aut <- ssprojection$p.aut
        projection$p.male.sex <- ssprojection$p.sex
    } else {
        fprojection <- .calcPoolProjection(lr.female.aut,lr.female.sex,opts)
        mprojection <- .calcPoolProjection(lr.male.aut,lr.male.sex,opts)
        projection$p.female.aut <- fprojection$p.aut
        projection$p.female.sex <- fprojection$p.sex
        projection$p.male.aut <- mprojection$p.aut
        projection$p.male.sex <- mprojection$p.sex
    }

    ## creating pool
    pool <- pd

    if (opts$pool.method %in% c("pca", "pca-nosex")) {
        pool$method <- "pca"
    } else if (opts$pool.method %in% c("ica", "ica-nosex")) {
        pool$method <- "ica"
    }

    pool$female$cov <- NULL
    pool$female$cov.med <- mc.female
    pool$female$projection.aut <- projection$p.female.aut
    pool$female$projection.sex <- projection$p.female.sex
    pool$female$filter <- filter.female
    pool$female$nvar <- nvar$female

    pool$male$cov <- NULL
    pool$male$cov.med <- mc.male
    pool$male$projection.aut <- projection$p.male.aut
    pool$male$projection.sex <- projection$p.male.sex
    pool$male$filter <- filter.male
    pool$male$nvar <- nvar$male

    return(pool)
}

.sexSpecificPool <- function(gpool, sex) {
    sspool <- gpool
    sspool$gender <- gpool[[sex]]
    sspool$male <- NULL
    sspool$female <- NULL
    sspool$sex <- sex
    sspool$projection.aut <- sspool$gender$projection.aut
    sspool$projection.sex <- sspool$gender$projection.sex
    sspool$cov.med <- sspool$gender$cov.med
    sspool$filter <- sspool$gender$filter
    sspool$gender <- NULL
    return(sspool)
}
