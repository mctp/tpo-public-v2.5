importClones <- function(clns.fns, fields="-nFeature CDR3 -count -fraction -targetQualities -vHit -vAlignment -jHit -jAlignment -cHit -cAlignment") {
    clone.tbls <- foreach(clns.fn=clns.fns) %do% {
            clone.tbl <- fread(paste(c(system2("mixcr", sprintf("exportClones %s -f %s /dev/stdout", fields, clns.fn),
                                    stdout=TRUE, stderr=FALSE), ""), collapse="\n"))
            clone.tbl$id <- basename(dirname(clns.fn))
            return(clone.tbl)
    }
    clone.tbl <- rbindlist(clone.tbls[sapply(clone.tbls, is.data.table)])
    return(clone.tbl)
}

makeClonotypes <- function(clone.tbl, vdjc) {
    tmp <- cbind(clone.tbl,
                   fread(paste(clone.tbl$bestVAlignment, collapse="\n"), sep="|", header=FALSE,
                         col.names=c("v.r.beg", "v.r.end", "v.r.len", "v.q.beg", "v.q.end", "mut", "score")),
                   fread(paste(clone.tbl$bestJAlignment, collapse="\n"), sep="|", header=FALSE,
                         col.names=c("j.r.beg", "j.r.end", "j.r.len", "j.q.beg", "j.q.end", "mut", "score")))
    tmp[,
    ":="(
        Vflank=str_sub(str_sub(vdjc[bestVHit]$seq, start=1, end=-(v.r.len-v.r.beg-19)), start=-40),
        Jflank=str_sub(vdjc[bestJHit]$seq, start=-(j.r.len-j.r.end))
    )
    ]
    tmp[,.(id, bestVHit, bestJHit, bestCHit, cloneCount, cloneFraction, Vflank, nSeqCDR3, Jflank)]
}

prepareVDJC <- function() {
    {
        fn <- system.file("extdata", "repseqio.v1.2.json", package="rnascape.api")
        db <- read_json(fn)[[1]] # Human
        hs.igx <- rbindlist(lapply(db[["genes"]], data.frame), fill=TRUE)
        hs.igx[,seqid:=str_match(baseSequence, "//(.*)")[,2]]
        hs.igx.fa <- readDNAStringSet("data/igx.fasta")
        names(hs.igx.fa) <- str_match(names(hs.igx.fa), "^[^\\s]+")[,1]
        hs.igx.ap <- rbind(
            hs.igx[geneType=="V", .(name, geneType, isFunctional, seqid, begin=anchorPoints.L2Begin+1, end=anchorPoints.VEnd)],
            hs.igx[geneType=="D", .(name, geneType, isFunctional, seqid, begin=anchorPoints.DBegin+1, end=anchorPoints.DEnd)],
            hs.igx[geneType=="J", .(name, geneType, isFunctional, seqid, begin=anchorPoints.JBegin+1, end=anchorPoints.FR4End)],
            hs.igx[geneType=="C", .(name, geneType, isFunctional, seqid, begin=anchorPoints.CBegin+1, end=anchorPoints.CExon1End)]
        )
        hs.igx.ap[,rc:=begin>end]
        hs.igx.ap.fwd <- hs.igx.ap[(!rc)]
        hs.igx.ap.fwd$seq <- as.character(subseq(hs.igx.fa[hs.igx.ap.fwd$seqid], hs.igx.ap.fwd$begin, hs.igx.ap.fwd$end))
        hs.igx.ap.rev <- hs.igx.ap[(rc)]
        hs.igx.ap.rev[,":="(begin=end+1L, end=begin-1L, rc=TRUE)]
        hs.igx.ap.rev$seq <- as.character(reverseComplement(subseq(hs.igx.fa[hs.igx.ap.rev$seqid],
                                                                   hs.igx.ap.rev$begin, hs.igx.ap.rev$end)))
        hs.igx.ap <- rbind(hs.igx.ap.fwd, hs.igx.ap.rev)
    }
    setkey(hs.igx.ap, name)
    return(hs.igx.ap)
}
