OPTS <- list(
    drop.chrM=TRUE,
    ## cin
    cin.AI.min.len = 2e6,
    cin.arm.level.portion = 0.8,
    cin.telomere.len = 1e4,
    cin.NtAI.min.len = 1e6,
    cin.LOH.min.len = 15e6,
    cin.lst.min.seglen = 10e6,
    cin.lst.min.segdist = 3e6,
    cin.SI.bucket.min = 0,
    cin.SI.bucket.max = 10,
    cin.SI.bucket.len = 0.5,
    cin.AiCNA.min.len = 8e6,
    cin.AbCNA.min.len = 8e6,
    cin.CnCNA.min.len = 4e6,
    ## ccf
    ccf.dp.af.max = 300
)
