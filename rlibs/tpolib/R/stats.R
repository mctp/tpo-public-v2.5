#' @export
row.fisher.test <- function(row, ...) {
    f <- fisher.test(matrix(row, nrow = 2), ...)
    return(c(
        p_val = f$p.value,
        or = f$estimate[[1]]
    ))
}

#' @export
log_sum_exp <- function(x) {
    offset <- max(x)
    log(sum(exp(x - offset))) + offset
}

#' @export
plog_sum_exp <- function(u, v) {
    m <- pmax(u, v)
    m + log(exp(u - m) + exp(v - m))
}

#' @export
max_na_rm <- function(x) max(x, na.rm=TRUE)

#' @export
absMedDiff <- function(x, y, na.rm=TRUE) {
    amd <- abs(median(x, na.rm=na.rm)-median(y, na.rm=na.rm))
    return(amd)
}

#' @export
minSum <- function(x, y, na.rm=TRUE) {
    ms <- min(sum(x, na.rm=na.rm), sum(y, na.rm=na.rm))
    return(ms)
}

#' @export
sigmoid <- function(x,xcord.top,xcord.low) {
    ## a function to move points onto a sigmoid function with middle to on zero
    xmiddle <- (xcord.top + xcord.low)/2
    xleft <- min(xcord.top,xcord.low)
    xright <- max(xcord.top,xcord.low)
    if (xcord.top < xcord.low) { #inverse
        sig <- 1/(1+exp(-(((2*xmiddle-x)-(xleft+xright)/2))/((xright - xleft)/2)*5))
    } else {
        sig <- 1/(1+exp(-((x-(xleft+xright)/2))/((xright - xleft)/2)*5))
    }
    return(sig)
}
