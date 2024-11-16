## hets filters
.filterGenomeNormalHetsSomatic <- function(var, common, opts) {
    mask.strict <- var$mask.strict | var$mask.giab
    lo_pass <- opts$baf.het.range.genome[1]
    hi_pass <- opts$baf.het.range.genome[2]
    ## genome germline hets
    hets <-
        ## germline SNV
        !var$SOMATIC &
        var$TYPE=="SNV" &
        ## heterozygous in normal
        var$n.GT %in% c("0/0/0/1", "0/0/1/1", "0/1/1/1", "0/1") &
        ## high quality
        var$n.DP > opts$baf.het.min.dp.genome &
        (
            (!mask.strict & !var$mask.loose &
             var$n.AF > lo_pass + 0.00 & var$n.AF < hi_pass - 0.00) |
            (mask.strict & !var$mask.loose &
             var$n.AF > lo_pass + 0.03 & var$n.AF < hi_pass - 0.03) |
            (mask.strict &  var$mask.loose &
             var$n.AF > lo_pass + 0.06 & var$n.AF < hi_pass - 0.06)
        ) &
        ## common variants
        common
    return(hets)
}

.filterGenomeNormalHetsGermline <- function(var, common, opts) {
    mask.strict <- var$mask.strict | var$mask.giab
    ## genome germline hets
    hets <-
        ## germline SNV
        !var$SOMATIC &
        var$TYPE=="SNV" &
        ## heterozygous in normal
        var$n.GT %in% c("0/0/0/1", "0/0/1/1", "0/1/1/1", "0/1") &
        ## high quality
        var$n.DP > opts$baf.het.min.dp.genome &
        !mask.strict &
        ## common variants
        common
    return(hets)
}

.filterTargetNormalHetsSomatic <- function(var, common, opts) {
    mask.strict <- var$mask.strict | var$mask.giab
    lo_pass <- opts$baf.het.range.target[1]
    hi_pass <- opts$baf.het.range.target[2]
    ## target germline hets
    hets <-
        ## germline SNV
        !var$SOMATIC &
        var$TYPE=="SNV" &
        ## heterozygous in normal
        var$n.GT %in% c("0/0/0/1", "0/0/1/1", "0/1/1/1", "0/1") &
        ## high quality
        var$n.DP > opts$baf.het.min.dp.target &
        (
            (!mask.strict & !var$mask.loose &
             var$n.AF > lo_pass + 0.00 & var$n.AF < hi_pass - 0.00) |
            (mask.strict & !var$mask.loose &
             var$n.AF > lo_pass + 0.03 & var$n.AF < hi_pass - 0.03) |
            (mask.strict &  var$mask.loose &
             var$n.AF > lo_pass + 0.06 & var$n.AF < hi_pass - 0.06)
        ) &
        ## common variants
        common
    return(hets)
}

.filterTargetNormalHetsGermline <- function(var, common, opts) {
    mask.strict <- var$mask.strict | var$mask.giab
    ## target germline hets
    hets <-
        ## germline SNV
        !var$SOMATIC &
        var$TYPE=="SNV" &
        ## heterozygous in normal
        var$n.GT %in% c("0/0/0/1", "0/0/1/1", "0/1/1/1", "0/1") &
        ## high quality
        var$n.DP > opts$baf.het.min.dp.target &
        !mask.strict &
        ## common variants
        common
    return(hets)
}

.filterGenomeTumorHets <- function(var, common, opts) {
    ## tumor-only heterozygous variants
    mask.strict <- var$mask.strict | var$mask.giab
    ## genome hets
    het <-
        ## germline SNV
        !var$SOMATIC &
        var$TYPE=="SNV" &
        ## high quality
        var$t.DP > opts$baf.het.min.dp.genome &
        !mask.strict &
        ## heterozygous in tumor
        (
          (var$t.GT %in% c("0/0/0/1", "0/0/1/1", "0/1/1/1", "0/1") & common) |
          (var$t.AF > 0.05 & var$t.AF < 0.95)
        )
    return(het)
}

.filterTargetTumorHets <- function(var, common, opts) {
    mask.strict <- var$mask.strict | var$mask.giab
    ## genome hets
    het <-
        ## germline SNV
        !var$SOMATIC &
        var$TYPE=="SNV" &
        ## high quality
        var$t.DP > opts$baf.het.min.dp.target &
        !mask.strict &
        ## heterozygous in tumor
        (
          (var$t.GT %in% c("0/0/0/1", "0/0/1/1", "0/1/1/1", "0/1") & common) |
          (var$t.AF > 0.05 & var$t.AF < 0.95)
        )
    return(het)
}

## Coverage filters
.filterGenomeTumorNormalCoverage <- function(var, opts) {
    t.hi.dp <- quantile(var$t.DP, 0.95, na.rm=TRUE)
    t.lo.dp <- quantile(var$t.DP, 0.05, na.rm=TRUE)
    n.hi.dp <- quantile(var$n.DP, 0.95, na.rm=TRUE)
    n.lo.dp <- quantile(var$n.DP, 0.05, na.rm=TRUE)
    ## coverage
    covered <-
        var$t.DP > t.lo.dp & var$t.DP < t.hi.dp &
        var$n.DP > n.lo.dp & var$n.DP < n.hi.dp &
        var$t.DP > opts$baf.cov.min.dp.genome
    return(covered)
}

.filterTargetTumorNormalCoverage <- function(var, opts) {
    t.hi.dp <- quantile(var$t.DP, 0.99, na.rm=TRUE)
    t.lo.dp <- quantile(var$t.DP, 0.01, na.rm=TRUE)
    n.hi.dp <- quantile(var$n.DP, 0.99, na.rm=TRUE)
    n.lo.dp <- quantile(var$n.DP, 0.01, na.rm=TRUE)
    ## coverage
    covered <-
        var$t.DP > t.lo.dp & var$t.DP < t.hi.dp &
        var$n.DP > n.lo.dp & var$n.DP < n.hi.dp &
        var$t.DP > opts$baf.cov.min.dp.target
    return(covered)
}

.filterGenomeTumorCoverage <- function(var, opts) {
    t.hi.dp <- quantile(var$t.DP, 0.95, na.rm=TRUE)
    t.lo.dp <- quantile(var$t.DP, 0.05, na.rm=TRUE)
    ## coverage
    covered <-
        var$t.DP > t.lo.dp & var$t.DP < t.hi.dp &
        var$t.DP > opts$baf.cov.min.dp.genome
    return(covered)
}

.filterTargetTumorCoverage <- function(var, opts) {
    t.hi.dp <- quantile(var$t.DP, 0.99, na.rm=TRUE)
    t.lo.dp <- quantile(var$t.DP, 0.01, na.rm=TRUE)
    ## coverage
    covered <-
        var$t.DP > t.lo.dp & var$t.DP < t.hi.dp &
        var$t.DP > opts$baf.cov.min.dp.target
    return(covered)
}

.filterGenomeNormalCoverage <- function(var, opts) {
    n.hi.dp <- quantile(var$n.DP, 0.95, na.rm=TRUE)
    n.lo.dp <- quantile(var$n.DP, 0.05, na.rm=TRUE)
    ## coverage
    covered <-
        var$n.DP > n.lo.dp & var$n.DP < n.hi.dp &
        var$n.DP > opts$baf.cov.min.dp.genome
    return(covered)
}

.filterTargetNormalCoverage <- function(var, opts) {
    n.hi.dp <- quantile(var$n.DP, 0.99, na.rm=TRUE)
    n.lo.dp <- quantile(var$n.DP, 0.01, na.rm=TRUE)
    ## coverage
    covered <-
        var$n.DP > n.lo.dp & var$n.DP < n.hi.dp &
        var$n.DP > opts$baf.cov.min.dp.target
    return(covered)
}

.filterGenomeTumorCommon <- function(var, opts) {
    common <- var$population.af > 0.025 & var$population.af < 0.975
    return(common)
}

.filterTargetTumorCommon <- function(var, opts) {
    common <- var$population.af > 0.01 & var$population.af < 0.99
    return(common)
}

.filterGenomeNormalCommon <- function(var, opts) {
    common <- var$population.af > 0.025 & var$population.af < 0.975
    return(common)
}

.filterTargetNormalCommon <- function(var, opts) {
    common <- var$population.af > 0.01 & var$population.af < 0.99
    return(common)
}


## On-target filters
.filterOnTarget <- function(tile, var, opts) {
    splash <- tile[tile$target] + opts$tile.shoulder
    ontarget <- var %over% splash
    return(ontarget)
}


## filters integration
.passTumorNormalHets <- function(tile, var, opts) {
    pass <- logical(length(var))
    if (any(var$SOURCE=="genome")) {
        var.genome <- var[var$SOURCE=="genome"]
        common.genome <- .filterGenomeNormalCommon(var.genome, opts)
        het.genome <- .filterGenomeNormalHetsSomatic(var.genome, common.genome, opts)
        covered.genome <- .filterGenomeNormalCoverage(var.genome, opts)
        pass[var$SOURCE=="genome"] <- het.genome & covered.genome
    }
    if (any(var$SOURCE=="target")) {
        var.target <- var[var$SOURCE=="target"]
        common.target <- .filterTargetNormalCommon(var.target, opts)
        het.target <- .filterTargetNormalHetsSomatic(var.target, common.target, opts)
        covered.target <- .filterTargetNormalCoverage(var.target, opts)
        ontarget <- .filterOnTarget(tile, var.target, opts)
        pass[var$SOURCE=="target"] <- het.target & covered.target & ontarget
    }
    pass[is.na(pass)] <- FALSE
    return(pass)
}

.passTumorOnlyHets <- function(tile, var, opts) {
    pass <- logical(length(var))
    if (any(var$SOURCE=="genome")) {
        var.genome <- var[var$SOURCE=="genome"]
        common.genome <- .filterGenomeTumorCommon(var.genome, opts)
        het.genome <- .filterGenomeTumorHets(var.genome, common.genome, opts)
        covered.genome <- .filterGenomeTumorCoverage(var.genome, opts)
        pass[var$SOURCE=="genome"] <- het.genome & covered.genome
    }
    if (any(var$SOURCE=="target")) {
        var.target <- var[var$SOURCE=="target"]
        common.target <- .filterTargetTumorCommon(var.target, opts)
        het.target <- .filterTargetTumorHets(var.target, common.target, opts)
        covered.target <- .filterTargetTumorCoverage(var.target, opts)
        ontarget <- .filterOnTarget(tile, var.target, opts)
        pass[var$SOURCE=="target"] <- het.target & covered.target & ontarget
    }
    pass[is.na(pass)] <- FALSE
    return(pass)
}

.passNormalOnlyHets <- function(tile, var, opts) {
    pass <- logical(length(var))
    if (any(var$SOURCE=="genome")) {
        var.genome <- var[var$SOURCE=="genome"]
        common.genome <- .filterGenomeNormalCommon(var.genome, opts)
        het.genome <- .filterGenomeNormalHetsGermline(var.genome, common.genome, opts)
        covered.genome <- .filterGenomeNormalCoverage(var.genome, opts)
        pass[var$SOURCE=="genome"] <- het.genome & covered.genome
    }
    if (any(var$SOURCE=="target")) {
        var.target <- var[var$SOURCE=="target"]
        common.target <- .filterTargetNormalCommon(var.target, opts)
        het.target <- .filterTargetNormalHetsGermline(var.target, common.target, opts)
        covered.target <- .filterTargetNormalCoverage(var.target, opts)
        ontarget <- .filterOnTarget(tile, var.target, opts)
        pass[var$SOURCE=="target"] <- het.target & covered.target & ontarget
    }
    pass[is.na(pass)] <- FALSE
    return(pass)
}

## filters API
passTumor <- function(tile, var, opts) {
    if ("n.GT" %in% names(mcols(var))) {
        pass <- .passTumorNormalHets(tile, var, opts)
    } else {
        pass <- .passTumorOnlyHets(tile, var, opts)
    }
    return(pass)
}

passNormal <- function(tile, var, opts) {
    pass <- .passNormalOnlyHets(tile, var, opts)
    return(pass)
}
