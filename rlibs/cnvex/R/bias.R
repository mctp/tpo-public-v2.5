.correctBias <- function(lr, target, hq, bias, biasdata, opts) {
  gt.essentials <- data.table(target = target, hq = hq, biasdata = biasdata)
  ## selecting the biasdata opts
  if (bias == "gc") { ## add these parameters to opts
      adjust.on <- opts$bias.gc.adjust.on
      adjust.off <- opts$bias.gc.adjust.off
      adjust.span.on <- opts$bias.gc.adjust.span.on
      adjust.span.off <- opts$bias.gc.adjust.span.off
      adjust.offset <- opts$bias.gc.adjust.offset
  ## } else if (bias == "meth") {
      ## TODO: BiasCorrection
      ## adjust.on <- opts$bias.meth.adjust.on
      ## adjust.off <- opts$bias.meth.adjust.off
      ## adjust.span.on <- opts$bias.meth.adjust.span.on
      ## adjust.span.off <- opts$bias.meth.adjust.span.off
      ## adjust.offset <- opts$bias.meth.adjust.offset
  ## } else if (bias == "reptime") {
      ## TODO: BiasCorrection
      ## adjust.on <- opts$bias.reptime.adjust.on
      ## adjust.off <- opts$bias.reptime.adjust.off
      ## adjust.span.on <- opts$bias.reptime.adjust.span.on
      ## adjust.span.off <- opts$bias.reptime.adjust.span.off
      ## adjust.offset <- opts$bias.reptime.adjust.offset
  }

  ## normalize, smooth, and bias-correct
  for (sel in unique(gt.essentials$target)) {
    lr.sel <- lr[gt.essentials$target==sel]
    gt.essentials.sel <- gt.essentials[gt.essentials$target==sel]
    weight.sel <- ifelse(gt.essentials.sel$hq, 1, 0)
    span.sel <- ifelse(sel, adjust.span.on, adjust.span.off)
    biasdata.residuals.sel <- limma::loessFit(y=lr.sel, x=gt.essentials.sel$biasdata,
                                              weight=weight.sel, min.weight=1e-9,
                                              span=span.sel)$residuals

    if (adjust.offset & !all(is.na(biasdata.residuals.sel))) {
      lr.offset.sel <- lm(biasdata.residuals.sel~1, weights=weight.sel)$coefficients[1]
      biasdata.residuals.sel <- biasdata.residuals.sel - lr.offset.sel
    }
    if (sel) {
      biasdata.range <- gt.essentials$biasdata > adjust.on[1] & gt.essentials$biasdata < adjust.on[2]
      biasdata.range[is.na(biasdata.range)] <- FALSE
      lr[gt.essentials$target & biasdata.range] <- biasdata.residuals.sel[biasdata.range[gt.essentials$target]]
    } else {
      biasdata.range <- gt.essentials$biasdata > adjust.off[1] & gt.essentials$biasdata < adjust.off[2]
      biasdata.range[is.na(biasdata.range)] <- FALSE
      lr[!gt.essentials$target & biasdata.range] <- biasdata.residuals.sel[biasdata.range[!gt.essentials$target]]
    }
  }
  lr[!is.finite(lr)] <- NA_real_
  return(lr)
}
