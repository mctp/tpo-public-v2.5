OPTS <- list(
    
    ## tiling
    tile.width=NA_integer_,
    tile.min.gap=50000,
    tile.shoulder=300,
    tile.hq.max.gap=0.005,
    tile.hq.min.unmasked=0.25,
    tile.hq.max.blacklist=0.001,
    tile.hq.min.totalcov=5000,
    tile.hq.max.giab.difficults=0.3,
    
    ## log-ratio
    lr.smooth="outlier",
    lr.smooth.window=13,
    lr.tumor='pair',
    lr.normal='pool',
    
    ## BAF
    baf.min.dp=6,
    baf.max.eff.dp=300,
    baf.het.range.genome=c(0.36, 0.64),
    baf.het.range.target=c(0.40, 0.60),
    baf.het.min.dp.genome=12,
    baf.het.min.dp.target=24,
    baf.cov.min.dp.genome=24,
    baf.cov.min.dp.target=48,
    
    ## segmentation
    seg.method="CBS",
    seg.only.target=TRUE,
    seg.cbs.baf=list(alpha=0.01, trim=0.025, min.width=2),
    seg.cbs.lr=list(alpha=0.01, trim=0.025, min.width=2),
    seg.cbs.weighted=FALSE,
    seg.rbs.selection="Lebarbier",
    
    ## pruning
    prune.lr.lo.threshold=1.5,
    prune.baf.lo.threshold=1.0,
    prune.len=TRUE,
    prune.lr.len.penalty=2.5,
    prune.baf.len.penalty=3.0,
    prune.nvar=TRUE,
    prune.lr.nvar.penalty=2.0,
    prune.baf.nvar.penalty=2.0,
    prune.hq=TRUE,
    prune.lr.hq.penalty=2.0,
    prune.baf.hq.penalty=2.0,
    
    ## Bias
    bias.gc.correct=TRUE,
    bias.gc.adjust.trend=TRUE,
    bias.gc.adjust.offset=TRUE,
    bias.gc.adjust.span.on=0.5,
    bias.gc.adjust.span.off=0.5,
    bias.gc.adjust.on=c(0.2, 0.8),
    bias.gc.adjust.off=c(NA_real_, NA_real_),
    ## bias.reptime.correct=FALSE,
    ## bias.reptime.adjust.trend=TRUE,
    ## bias.reptime.adjust.offset=TRUE,
    ## bias.reptime.adjust.span.on=0.5,
    ## bias.reptime.adjust.span.off=0.5,
    ## bias.reptime.adjust.on=c(0, 1),
    ## bias.reptime.adjust.off=c(NA_real_, NA_real_),
    ## bias.meth.correct=FALSE,
    ## bias.meth.adjust.trend=TRUE,
    ## bias.meth.adjust.offset=TRUE,
    ## bias.meth.adjust.span.on=0.5,
    ## bias.meth.adjust.span.off=0.5,
    ## bias.meth.adjust.on=c(0, 1),
    ## bias.meth.adjust.off=c(NA_real_, NA_real_),
    
    ## pool
    pool.method="pca",
    pool.lo.cov=0.25,
    pool.hi.cov=4,
    pool.hi.zero=0.20,
    pool.hi.nvar=0.95,
    pool.lo.nvar=NA_real_,
    pool.sd.out=3,
    pool.n.comp=20,
    pool.filter=FALSE,        # TRUE:filter out noise-prone tiles | FALSE: use all targeted tiles
    
    ## llik
    opt.max.C=24,
    opt.max.sC=36,
    opt.max.len.per.probe=1e6,
    opt.max.snp.per.segment=Inf,
    opt.p.af.anom=0.001,
    opt.p.lr.anom=0.01,
    opt.dp.af.max=300,
    opt.grid.llik="heuristic-v2",
    opt.grid.subc.exp=1,
    opt.grid.subc.prob=0.5,
    opt.only.hq=TRUE,
    
    ## grid
    opt.p.lo=0.10,
    opt.p.hi=0.99,
    opt.P.lo=1,
    opt.P.hi=6,
    opt.grid.max.C=12,
    opt.grid.p.res=0.0125,
    opt.fine.p.res=0.0050,
    opt.fine.p.off=0.0250,
    opt.cand.res=0.075,
    opt.cand.max.iter=5,
    
    ## sex detection
    sex.x.snpratio.threshold=0.25,
    sex.x.covratio.threshold=0.75,
    sex.y.covzero.threshold=0.05,
    
    ## output
    output.col.name=FALSE,
    output.same.marker=TRUE,

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
    cin.CnCNA.min.len = 4e6
)
