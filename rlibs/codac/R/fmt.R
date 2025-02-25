#' @export
bsFormat <- function(tbl, reorder=TRUE) {
    rep <- tbl[,.(
        gene.5=ifelse(gene_names.5.1!="",                 str_match(gene_names.5.1, "[^:]*")[,1],
               ifelse(gene_names.3.2!="", sprintf("(%s)", str_match(gene_names.3.2, "[^:]*")[,1]), cyt.5.1)),
        junction.5=sprintf("%s:%s:%s", chr.5, pos.5, str.5),
        gene.3=ifelse(gene_names.3.1!="",                 str_match(gene_names.3.1, "[^:]*")[,1],
               ifelse(gene_names.5.2!="", sprintf("(%s)", str_match(gene_names.5.2, "[^:]*")[,1]), cyt.3.1)),
        junction.3=sprintf("%s:%s:%s", chr.3, pos.3, str.3),
        spanning.reads=sum.jnc,
        breakpoints=sum.bpt,
        distance=tx.dst,
        spliced=d2a,
        "cpm 5'"=format(round(cpm.5, 2), nsmall=2, scientific=999),
        "cpm 3'"=format(round(cpm.3, 2), nsmall=2, scientific=999),
        "read fraction 5'" = format(round(sum.jnc / (tot.jnc.5 + tot.sp.jnc.5), 2), nsmall=2, scientific=999),
        "read fraction 3'" = format(round(sum.jnc / (tot.jnc.3 + tot.sp.jnc.3), 2), nsmall=2, scientific=999),
        "repetitive(5';3)'" = paste(art.5, art.3, sep=";"),
        chain=bs.chain
    )]
    if (reorder) {
        rep <- rep[order(chain)]
    }
    return(rep)
}

#' @export
slFormat <- function(tbl, reorder=TRUE) {
    rep <- tbl[,.(
        gene.5=ifelse(gene_names.5.1!="",                 str_match(gene_names.5.1, "[^:]*")[,1],
               ifelse(gene_names.3.2!="", sprintf("(%s)", str_match(gene_names.3.2, "[^:]*")[,1]), cyt.5.1)),
        junction.5=sprintf("%s:%s:%s", chr.5, pos.5, str.5),
        gene.3=ifelse(gene_names.3.1!="",                 str_match(gene_names.3.1, "[^:]*")[,1],
               ifelse(gene_names.5.2!="", sprintf("(%s)", str_match(gene_names.5.2, "[^:]*")[,1]), cyt.3.1)),
        junction.3=sprintf("%s:%s:%s", chr.3, pos.3, str.3),
        spanning.reads=sum.jnc,
        breakpoints=sum.bpt,
        spliced=ifelse(d2a, "✓", ""),
        HQ=ifelse(hq.bpt, "✓", ""),
        inframe=ifelse(orf, "✓", ""),
        distance=dst,
        topology=topo,
        "cpm 5'"=format(round(cpm.5, 2), nsmall=2, scientific=999),
        "cpm 3'"=format(round(cpm.3, 2), nsmall=2, scientific=999),
        "repetitive(5';3')" = paste(art.5, art.3, sep=";"),
        chain=sl.chain
    )]
    if (reorder) {
        rep <- rep[order(chain)]
    }
    return(rep)
}

#' @export
svFormat <- function(tbl, reorder=TRUE) {
    rep <- tbl[,.(
        gene.5=ifelse(gene_names.5.1!="",                 str_match(gene_names.5.1, "[^:]*")[,1],
               ifelse(gene_names.3.2!="", sprintf("(%s)", str_match(gene_names.3.2, "[^:]*")[,1]), cyt.5.1)),
        junction.5=sprintf("%s:%s:%s", chr.5, pos.5, str.5),
        gene.3=ifelse(gene_names.3.1!="",                 str_match(gene_names.3.1, "[^:]*")[,1],
               ifelse(gene_names.5.2!="", sprintf("(%s)", str_match(gene_names.5.2, "[^:]*")[,1]), cyt.3.1)),
        junction.3=sprintf("%s:%s:%s", chr.3, pos.3, str.3),
        spanning.reads=sum.jnc * (l1=="spn"),
        encompassing.reads=tot.enc,
        total.reads=tot.jnc+tot.enc,
        breakpoints=sum.bpt,
        spliced=ifelse(d2a, "✓", ""),
        HQ=ifelse(hq.bpt, "✓", ""),
        TS=ifelse(ts.warn, "✓", ""),
        inframe=ifelse(orf, "✓", ""),
        distance=dst,
        topology=topo,
        mm2.valid=paste(ifelse(mm2.valid, "✓", "")),
        gmap.valid=paste(ifelse(gmap.valid, "✓", "")),
        "read fraction 5'" = format(round(sum.jnc / (tot.jnc.5 + tot.sp.jnc.5), 2), nsmall=2, scientific=999),
        "read fraction 3'" = format(round(sum.jnc / (tot.jnc.3 + tot.sp.jnc.3), 2), nsmall=2, scientific=999),
        "recurrent(5';3')" = paste(unq.rec.5, unq.rec.3, sep=";"),
        "repetitive(5';3')" = paste(art.5, art.3, sep=";"),
        chain=sv.chain,
        contig=sapply(lapply(ctg.seq, as.character), "[", 1)
    )]
    if (reorder) {
        rep <- rep[order(chain)]
    }
    return(rep)
}

#' @export
tsFormat <- function(tbl, reorder=TRUE) {
    rep <- tbl[,.(
        gene.5=ifelse(gene_names.5.1!="",                 str_match(gene_names.5.1, "[^:]*")[,1],
               ifelse(gene_names.3.2!="", sprintf("(%s)", str_match(gene_names.3.2, "[^:]*")[,1]), cyt.5.1)),
        junction.5=sprintf("%s:%s:%s", chr.5, pos.5, str.5),
        gene.3=ifelse(gene_names.3.1!="",                 str_match(gene_names.3.1, "[^:]*")[,1],
               ifelse(gene_names.5.2!="", sprintf("(%s)", str_match(gene_names.5.2, "[^:]*")[,1]), cyt.3.1)),
        junction.3=sprintf("%s:%s:%s", chr.3, pos.3, str.3),
        spanning.reads=sum.jnc,
        breakpoints=sum.bpt,
        spliced=ifelse(d2a, "✓", ""),
        HQ=ifelse(hq.bpt, "✓", ""),
        inframe=ifelse(orf, "✓", ""),
        distance=dst,
        topology=topo,
        "read fraction 5'" = format(round(sum.jnc / (tot.jnc.5 + tot.sp.jnc.5), 2), nsmall=2, scientific=999),
        "read fraction 3'" = format(round(sum.jnc / (tot.jnc.3 + tot.sp.jnc.3), 2), nsmall=2, scientific=999),
        "recurrent(5';3')" = paste(unq.rec.5, unq.rec.3, sep=";"),
        "repetitive(5';3')" = paste(art.5, art.3, sep=";"),
        chain=ts.chain
    )]
    if (reorder) {
        rep <- rep[order(chain)]
    }
    return(rep)
}

#' @export
neoFormat <- function(tbl) {
    rep <- tbl[,.(
        gene.5=ifelse(!is.na(gene_name.5), gene_name.5, "-"),
        gene.3=ifelse(!is.na(gene_name.5), gene_name.3, "-"),
        rna_seq=sapply(lapply(x.rna.3, as.character), "[", 1),
        start_codon=cds1,
        protein_seq=sapply(lapply(x.prot.3, as.character), "[", 1),
        start_chimeric=x.prot.pos.3,
        contig_len=ctg.len,
        contig_cov=ctg.cov,
        contig_seq=sapply(lapply(ctg.seq, as.character), "[", 1),
        inframe=inframe.3,
        homology=homology
    )]
    return(rep)
}
