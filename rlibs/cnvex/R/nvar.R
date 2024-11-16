.normalCovCorrectGC <- function(normal.cov, gc) {
    normal.cov.gcCorrect <- foreach(patient = 1:nrow(normal.cov)) %dopar% {
        gc.residuals <- limma::loessFit(y = normal.cov[patient,], x=gc)$residuals
        cov.offset <- lm(gc.residuals~normal.cov[patient,])$coefficients[1]
        gc.residuals <- gc.residuals - cov.offset
        return(gc.residuals)
    }
    return(t(do.call(cbind,normal.cov.gcCorrect)))
}

.calcVariance <- function(normal.cov, opts) {
    nvar <- sapply(X = 1:ncol(normal.cov), function(X, NCOV = normal.cov) {
        return(sd(NCOV[,X], na.rm = TRUE))
    })
    return(nvar)
}

.calcNvar <- function(pd, opts) {
    ## female
    normal.cov.female <- t(pd$female$cov)
    normal.variance.female <- .calcVariance(normal.cov.female, opts)

    ## male
    normal.cov.male <- t(pd$male$cov)
    normal.variance.male <- .calcVariance(normal.cov.male, opts)

    normal.variance <- list(male = as.numeric(normal.variance.male),female = as.numeric(normal.variance.female))
    return(normal.variance)
}

.addNvarToPd <- function(pd, nvar, opts) {
    pd$male$nvar <- nvar$male
    pd$female$nvar <- nvar$female
    return(pd)
}
