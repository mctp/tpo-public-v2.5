.variantRepr <- function(tmp) {
    tbl <- tmp[,.(
        var_id=ID,
        chr, pos, ref, alt,
        ## support
        aft=AFT, adt=ADT, dpt=DPT, afn=AFN, adn=ADN, dpn=DPN,
        ## quality
        tlod, nlod, xfet, str,
        ## genes
        gene_id=GENE, transcript_id=TRANSCRIPT,
        ## consequence
        hgvsp=HGVSp, hgvsc=HGVSc, consequence=Consequence,
        cosmic_cnt, gnomad_af, clinvar
    )]
    return(tbl)
}

somaticRepr <- function(tbl) {
    .variantRepr(tbl)
}

germlineRepr <- function(tbl) {
    .variantRepr(tbl)
}

structuralRepr <- function(tbl) {
    repr <- tbl[,.(
        var_id=bnd.id,
        chr1,pos1,str1=STRAND1,chr2,pos2,str2=STRAND2,
        ## support
        aft=AFT,adt=ADT,dpt=DPT,
        ## quality
        qual=QUAL,
        score=MANTA.SCORE,
        filter=MANTA.FILTER,
        ## genes
        gene_id1=GENE1,gene_id2=GENE2,
        ## consequences
        topology=topo,
        csq1=CSQ1,csq2=CSQ2,
        exons1=EXONS1,exons2=EXONS2,
        ## interpretation
        insert=insert,
        contig=MANTA.CONTIG
    )]
    return(repr)
}

fusionRepr <- function(tmp) {
    repr <- tmp[,.(
        var_id,
        chr5=chr.5, pos5=pos.5, str5=str.5, chr3=chr.3, pos3=pos.3, str3=str.3,
        ## support
        spn_reads=sum.jnc * (l1=="spn"),
        enc_reads=tot.enc,
            ## quality
        hq=hq.bpt,
        mm2=mm2.valid,
        gmap=gmap.valid,
        ## genes
        gene_ids.5.1, gene_ids.3.1,
        gene_ids.5.2, gene_ids.3.2,
        ## consequence
        distance=dst,
        topology=topo,
        spliced=d2a,
        inframe=orf,
        ## interpretation
        chain=sv.chain,
        contig=sapply(lapply(ctg.seq, as.character), "[", 1)
    )]
    return(repr)
}
