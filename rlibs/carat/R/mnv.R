.toGappedGR <- function(recs, bsgenome) {
    ## constructs a genomic range object with
    ## REF and ALT bases, including gaps between variants
    gr <- unname(granges(recs)[,c("REF", "ALT")])
    ## keep only real gaps
    gp <- gaps(gr)
    gp <- gp[seqnames(gp) == unique(seqnames(gr)) & strand(gp)==unique(strand(gr))]
    gp <- tail(head(gp, -1), -1)
    ##
    if (length(gp)>0) {
        mcols(gp)$REF <- getSeq(bsgenome, gp)
        mcols(gp)$ALT <- DNAStringSetList(as.list(as.character(mcols(gp)$REF)))
        mnv.gr <- sort(c(gr, gp))
    } else {
        mnv.gr <- gr
    }
    return(mnv.gr)
}

.findClumps <- function(recs, maxgap) {
    ## finds connected components of variants at most maxgap bases apart
    var.overlaps <- findOverlaps(recs, maxgap=maxgap)
    ## igraph
    clump.id <- clusters(graph_from_edgelist(as.matrix(var.overlaps)))$membership
    var.clumps <- split(recs, clump.id)
    return(var.clumps)
}

.toMNV <- function(recs, bsgenome) {
    ## convert a set of phased variants into a single block MNV
    mnv.gr <- .toGappedGR(recs, bsgenome)
    mnv <- sort(recs)[1] # select first variant has correct position
    ## correct REF,ALT,width of the mnv
    fixed(mnv)$REF <- DNAStringSet(paste0(mnv.gr$REF, collapse=""))
    fixed(mnv)$ALT <- DNAStringSetList(paste0(unlist(mnv.gr$ALT), collapse=""))
    width(mnv) <- width(ref(mnv))
    names(mnv) <- sprintf("%s:%s_%s/%s", seqnames(mnv), start(mnv), ref(mnv), unlist(alt(mnv)))
    return(mnv)
}

#' @export
phaseMNVSingle <- function(var, maxgap, bsgenome) {
    ## partition variants into phased and unphased
    is.phased <- !is.na(geno(var)$PID[,1])
    var.unphased <- var[!is.phased]
    var.phased <- var[is.phased]
    ## fix phase-ID (PID) to include chromosome i.e. to be genome-wide unique
    ## is it possible that PID[,2] was different from PID[,1]
    geno(var.phased)$PID[,1] <- paste(seqnames(var.phased), geno(var.phased)$PID[,1], sep = "_")
    ## split variants into phased groups (list of VCF objects)
    var.phased.split <- split(var.phased, geno(var.phased)$PID[,1])
    ## further split phased variants to no more than maxgap bases apart
    var.clumps.split <- unlist(List(lapply(var.phased.split, .findClumps,
                                             maxgap=maxgap)))
    ## select singleton clumps (i.e. not MNVs)
    var.phased.distal <- unlist(unname(var.clumps.split[lengths(var.clumps.split)==1]))
    ## select MNVs
    tmp <- var.clumps.split[lengths(var.clumps.split)!=1]
    ## convert block of phases SNVs/Indels into a single MNV

    var.phased.mnv <- unlist(List(unname(lapply(tmp, .toMNV,
                                                  bsgenome=bsgenome))))
    ## this is a not a hack anymore, rbind on VCF objects does not work
    vars <- list(var.unphased, var.phased.distal, var.phased.mnv)

    ## write vcfs to disk
    vcfs <- sapply(vars, function(var) {
        if(!is.null(var)) {
            out <- tempfile()
            out.gz <- paste0(out, ".gz")
            writeVcf(sort(var), out)
            system2("bgzip", out)
            system2("tabix", out.gz)
            return(out.gz)
        }
    })

    vcfs <- vcfs[!sapply(vcfs,is.null)]

    return(vcfs)
}

#' @export
phaseMNV <- function(var, maxgap, gobj) {
    ## split variants by chromosome
    var.split.chr <- split(var, seqnames(var), drop=FALSE)

    ## run phaseMNV function on each resulting vcf
    var.mnv.list <- foreach(x=var.split.chr) %dopar% {
        phaseMNVSingle(x, maxgap, gobj$seq)
    }

    ## combine into a single list
    vars <- unlist(var.mnv.list)

    return(vars)
}
