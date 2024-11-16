#' @export
merge_list <- function(x, y) {
  if (length(x) == 0)
    return(y)
  if (length(y) == 0)
    return(x)
  i <- match(names(y), names(x))
  i <- is.na(i)
  if (any(i))
    x[names(y)[which(i)]] <- y[which(i)]
  return(x)
}

#' @export
dtRound <- function(dt) {
  for (i in names(dt)) {
    if (class(dt[,get(i)])=="numeric") {
      dt[,":="((i), round(get(i), 6))]
    }
  }
  return(dt)
}

#' @export
naTo0 <- function(DT, omit=NULL, zero=0) {
  for (i in setdiff(names(DT), omit))
    DT[is.na(get(i)),(i):=zero]
  return(DT)
}

#' @export
naToBlank <- function(DT, omit=NULL) {
  for (i in setdiff(names(DT), omit))
    DT[is.na(get(i)),(i):=""]
  return(DT)
}

#' @export
omitBlank <- function(x) {
  return(x[!(x=="")])
}

#' @export
pickExisting <- function(...) {
  args <- list(...)
  for (arg in args) {
    if (!is.null(arg) && !is.na(arg)) {
      return(arg)
    }
  }
  return(NA)
}
