.resolveIds <- function(csq) {
    csq[,id.short:=NA_character_]
    csq[Consequence=="intergenic_variant",id.short:=paste0("intergenic_", ID)]
    csq[SOURCE=="RefSeq",id.short:=str_match(Feature, "(.*)\\.[0-9]+$")[,2]]
    csq[SOURCE=="Ensembl" | BIOTYPE=="aligned_transcript",id.short:=Feature]
    csq[,HGVSc.short:=str_match(HGVSc, "[^:]*$")[,1]]
    csq[,HGVSp.short:=str_match(HGVSp, "[^:]*$")[,1]]
    csq[,HGVS_Equivalent:=paste(id.short, collapse=";"),.(ID, HGVSc.short, HGVSp.short)]
    return(csq)
}

.impactScore <- function(csq) {
    csq$impact.score <- 0
    csq[IMPACT=="LOW", impact.score:=10]
    csq[IMPACT=="MODERATE", impact.score:=20]
    csq[IMPACT=="HIGH", impact.score:=30]
    csq[,impact.score:=impact.score +
            ifelse(EXON!="" & Protein_position!="" & BIOTYPE=="protein_coding", 7,
            ifelse(EXON!="" & cDNA_position!="" & BIOTYPE=="protein_coding", 5,
            ifelse(INTRON!="" & BIOTYPE=="protein_coding", 3, 0)))
    ]
    return(csq)
}

.preferenceScore <- function(csq, transcripts) {
    if (transcripts=="onco1500-v6a") {
        ## update w/ transcript preference data
        v6a.pref <- fread(system.file("extdata/transcripts/onco1500-v6a/transcript_preference.csv", package="carat"))
        setkey(v6a.pref, transcript_id)
        ## ranks to scores
        csq[,preference.score:=0]
        setkey(csq, id.short)
        csq[v6a.pref, ":="(Gene=ensg, rank=rank)]
        setkey(csq, ID)
        csq[rank==1, preference.score:=60]
        csq[rank>1, preference.score:=30]
    } else {
        ## update w/ transcript preference data
        known.transcripts <- fread(system.file("extdata/transcripts/generic/transcript_preference.csv", package="carat"))
        setkey(csq, id.short)
        setkey(known.transcripts, id.short)
        csq <- known.transcripts[csq]
        for (i in setdiff(colnames(known.transcripts), "id.short"))
            csq[is.na(get(i)),(i):=FALSE]
        setkey(csq, ID)
        ## make sure does not fail if columns are missing
        for (col in c("MANE_SELECT", "CANONICAL", "PICK", "CCDS", "TSL")) {
            if (!(col %in% colnames(csq))) {
                csq[,(col):=""]
            }
        }
        ## score preference of transcript
        csq[,preference.score:=0]
        csq[(known), preference.score:=3]
        csq[(known & BIOTYPE=="protein_coding" & FLAGS==""), preference.score:=5]
        ## RefSeq
        csq[SOURCE=="RefSeq" & mo.v4, preference.score:=10]
        csq[SOURCE=="RefSeq" & (ccds | lrg), preference.score:=20]
        csq[SOURCE=="RefSeq" & CANONICAL=="YES" & (mo.v4 | ccds | lrg), preference.score:=30]
        csq[SOURCE=="RefSeq" & select, preference.score:=40]
        csq[SOURCE=="RefSeq" & preference.score >=30 & PICK=="1", preference.score:=50]
        csq[SOURCE=="RefSeq" & preference.score >=20 & mo.v6, preference.score:=60]
        csq[SOURCE=="RefSeq" & MANE_SELECT!="", preference.score:=70]
        ## ENSEMBL
        csq[, valid.symbol:=!(SYMBOL_SOURCE %in% "Clone_based_ensembl_gene")]
        csq[SOURCE=="Ensembl" & valid.symbol & (basic | mo.v4), preference.score:=10]
        csq[SOURCE=="Ensembl" & valid.symbol & (lrg | ccds | CCDS!=""), preference.score:=20]
        csq[SOURCE=="Ensembl" & valid.symbol & CANONICAL=="YES" & (basic | mo.v4 | lrg | ccds | CCDS!=""), preference.score:=30]
        csq[SOURCE=="Ensembl" & valid.symbol & (appris1 | tsl1 | TSL=="1" | valid), preference.score:=40]
        csq[SOURCE=="Ensembl" & valid.symbol & preference.score >=30 & PICK=="1", preference.score:=50]
        csq[SOURCE=="Ensembl" & valid.symbol & preference.score >=20 & mo.v6, preference.score:=60]
        csq[SOURCE=="Ensembl" & MANE_SELECT!="", preference.score:=70]
        csq[SOURCE=="Ensembl", preference.score:=preference.score+1]
    }
    ## TODO: Warnings
    csq[, Warnings:=""]
    return(csq)
}

.scoreCSQ <- function(csq, transcripts) {
    csq <- .resolveIds(csq)
    csq <- .impactScore(csq)
    csq <- .preferenceScore(csq, transcripts)
    return(csq)
}

.pickCSQ2 <- function(var, bnd.pair, transcripts) {
    ## collect consequences (CSQ) for each variant breakend (BND)
    csq <- .csqTable(var)
    ## score each CSQ for each BND
    csq.score.full <- .scoreCSQ(csq, transcripts)
    csq.score <- csq.score.full[,.(ID, SYMBOL, id.short, impact.score, preference.score)]
    setkey(csq.score, ID)
    ## create table with scores for both breakends and all transcripts
    setkey(bnd.pair, id1)
    bnd.score.1 <- csq.score[bnd.pair]
    setkey(bnd.pair, id2)
    bnd.score.2 <- csq.score[bnd.pair]
    bnd.score.1 <- bnd.score.1[,.(bnd.id, proximal, SYMBOL1=SYMBOL, id.short, impact.score.1=impact.score, preference.score.1=preference.score)]
    bnd.score.2 <- bnd.score.2[,.(bnd.id, proximal, SYMBOL2=SYMBOL, id.short, impact.score.2=impact.score, preference.score.2=preference.score)]
    setkey(bnd.score.1, bnd.id, proximal, id.short)
    setkey(bnd.score.2, bnd.id, proximal, id.short)
    bnd.score <- merge(bnd.score.1, bnd.score.2, all=TRUE)
    ##
    setkey(bnd.score, bnd.id)
    setkey(bnd.pair, bnd.id)
    bnd.score[bnd.pair, ":="(id1=id1,id2=id2)]
    bnd.score[,pair.impact.score:=NA_integer_]
    bnd.score[,pair.preference.score:=NA_integer_]
    bnd.score[SYMBOL1=="", SYMBOL1:=NA_character_]
    bnd.score[,ISNA1:=all(is.na(SYMBOL1[preference.score.1>=3])),bnd.id]
    bnd.score[SYMBOL2=="", SYMBOL2:=NA_character_]
    bnd.score[,ISNA2:=all(is.na(SYMBOL2[preference.score.1>=3])),bnd.id]
    bnd.score[(SYMBOL1==SYMBOL2) | (proximal & (ISNA1 | ISNA2)), pair.impact.score:=pmax(impact.score.1, impact.score.2, na.rm=TRUE)]
    bnd.score[(SYMBOL1==SYMBOL2) | (proximal & (ISNA1 | ISNA2)), pair.preference.score:=pmax(preference.score.1, preference.score.2, na.rm=TRUE)]
    bnd.score.1 <- bnd.score[,.(ID=id1, id.short, pair.impact.score, pair.preference.score)]
    bnd.score.2 <- bnd.score[,.(ID=id2, id.short, pair.impact.score, pair.preference.score)]
    bnd.score.12 <- rbind(bnd.score.1, bnd.score.2)
    setkey(csq.score.full, ID, id.short)
    setkey(bnd.score.12, ID, id.short)
    csq.score.expand <- csq.score.full[bnd.score.12]
    ## pick
    csq.pick <- csq.score.expand[order(pair.preference.score<3, -pair.impact.score, -pair.preference.score,
                                       preference.score<3, -impact.score, -preference.score
                                       ),.SD[1],.(ID)]
    setkey(csq.pick, ID)
    return(csq.pick)
}

.pickCSQ1 <- function(var, transcripts) {
    ## collect consequences (CSQ) for each variant
    csq <- .csqTable(var)
    ## score consequences
    csq.score <- .scoreCSQ(csq, transcripts)
    ## Pick
    csq.pick <- csq.score[order(preference.score<3, -impact.score, -preference.score),.SD[1],.(ID)]
    setkey(csq.pick, ID)
    csq.pick <- csq.pick[names(var)]
    return(csq.pick)
}

.varRNG <- function(var) {
    rng <- granges(var)
    var.rng <- data.table(ID=names(rng), chr=as.character(seqnames(rng)), pos=start(rng),
                          ref=as.character(rng$REF), alt=as.character(rng$ALT), filter=as.character(rng$FILTER))
    var.rng$context <- info(var)$context
    var.rng$sbs96 <- info(var)$SBS96
    return(var.rng)
}

.varCSQ <- function(var, transcripts) {
    csq.pick <- .pickCSQ1(var, transcripts)
    csq.pick.short <- csq.pick[,.(
        ID,
        SYMBOL, GENE=Gene, TRANSCRIPT=id.short, EXON,
        HGVSc=HGVSc.short,
        HGVSp=HGVSp.short,
        Consequence, IMPACT, SIFT, Warnings, HGVS_Equivalent
    )]
    return(csq.pick.short)
}

.varANN <- function(var) {
    ID <- names(var)

    ## variant annotation
    dbsnp <- info(var)$dbSNP_RS
    clinvar <- sapply(info(var)$CLINVAR_CLNSIG, "[", 1)
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinid}/
    clinid <- sapply(info(var)$CLINVAR_ALLELEID, "[", 1)

    if (!is.null(info(var)$COSMIC_CNT)) {
        ## this is an IntegerList because COSMIC_CNT is an A
        stopifnot(all(lengths(info(var)$COSMIC_CNT)==1))
        cosmic_cnt <- unlist(info(var)$COSMIC_CNT)
        cosmic_cnt[is.na(cosmic_cnt)] <- 0L
    } else {
        cosmic_cnt <- NA_integer_
    }
    if (!is.null(info(var)$GNOMAD_AF)) {
        ## this is a NumericList because GNOMAD_AF is an A
        stopifnot(all(lengths(info(var)$GNOMAD_AF)==1))
        gnomad_af <- unlist(info(var)$GNOMAD_AF)
        gnomad_af[is.na(gnomad_af)] <- 0
    } else {
        gnomad_af <- NA_real_
    }
    ## this a vector because KG_AF is a 1
    if (!is.null(info(var)$KG_AF)) {
        kg_af <- info(var)$KG_AF
        kg_af[is.na(kg_af)] <- 0
    } else {
        kg_af <- NA_real_
    }
    ## this is a CharacterList becase PON is a .
    if (!is.null(info(var)$PON)) {
        pon <- as.list(info(var)$PON)
    } else {
        pon <- rep(list(character(0)), length(var))
    }
    ## this is a CharacterList with character(0) as it is a .
    ## rarely it can have multiple repeats not per allele
    rmsk_hit <- as.list(info(var)$RMSK_HIT)

    var.ann <- data.table(
        ID, dbsnp, clinid, clinvar, cosmic_cnt, gnomad_af, kg_af, pon, rmsk_hit
    )
    return(var.ann)
}

.varQNT <- function(var) {
    ID <- names(var)

    ## variant reads support
    AFT <- geno(var)$AF[,1]
    ADT <- geno(var)$AD[,1,2]
    DPT <- geno(var)$DP[,1]
    if (ncol(geno(var)$DP)==2) {
        AFN <- geno(var)$AF[,2]
        ADN <- geno(var)$AD[,2,2]
        DPN <- geno(var)$DP[,2]
    } else {
        AFN <- NA_integer_
        ADN <- NA_integer_
        DPN <- NA_integer_
    }
    if (is.null(geno(var)$ALT_F1R2)) {
        adt_fwd <- NA_integer_
    } else {
        adt_fwd <- geno(var)$ALT_F1R2[,1]
    }
    if (is.null(geno(var)$ALT_F2R1)) {
        adt_rev <- NA_integer_
    } else {
        adt_rev <- geno(var)$ALT_F2R1[,1]
    }
    if (is.null(geno(var)$QSS)) {
        mqs <- NA_real_
    } else {
        if (is(geno(var)$QSS, "matrix")) {
            mqs <- do.call(rbind, geno(var)$QSS[,1])[,2] / ADT
        } else { # sentieon version >= 202010.01
            mqs <- geno(var)$QSS[,1,2] / ADT
        }
    }

    ## variant probability
    if (is.null(info(var)$TLOD)) {
        tlod <- NA_real_
    } else {
        tlod <- unlist(info(var)$TLOD)
    }
    if (is.null(info(var)$NLOD)) {
        nlod <- NA_real_
    } else {
        nlod <- unlist(info(var)$NLOD)
    }
    if (is.null(info(var)$NLODF)) {
        flod <- NA_real_
    } else {
        flod <- info(var)$NLODF
    }
    if (is.null(info(var)$PV)) {
        xfet <- NA_real_
    } else {
        xfet <- info(var)$PV
    }
    if (is.null(info(var)$SOR)) {
        sor <- NA_real_
    } else {
        sor <- info(var)$SOR
    }
    if (is.null(geno(var)$MQRankSumPS)) {
        mqrs <- NA_real_
    } else {
        mqrs <- geno(var)$MQRankSumPS[,1]
    }
    if (is.null(info(var)$STR)) {
        str <- NA_integer_
    } else {
        str <- info(var)$STR
    }
    if (is.null(info(var)$RU) || is.null(info(var)$RPA)) {
        str_ru <- NA_character_
        str_len <- NA_integer_
        str_diff <- NA_integer_
    } else {
        str_ru <- info(var)$RU
        rpa <- as.list(info(var)$RPA)
        rpa[lengths(rpa)==0] <- list(c(0,0))
        rpa <- do.call(rbind, rpa)
        str_len <- rpa[,1]
        str_diff <- rpa[,1] - rpa[,2]
    }
    if (is.null(geno(var)$FOXOG)) {
        oxog <- NA_real_
    } else {
        oxog <- geno(var)$FOXOG[,1]
    }
    if (is.null(info(var)$ECNT)) {
        ecnt <- NA_real_
    } else {
        ecnt <- info(var)$ECNT
    }
    if (is.null(info(var)$CPTAC3_ADX)) {
        c3_adx <- NA_real_
    } else {
        c3_adx <- unlist(info(var)$CPTAC3_ADX)
    }
    if (is.null(geno(var)$PID)) {
        pid <- NA_real_
    } else {
        pid <- geno(var)$PID[,1]
    }

    ## output
    var.qnt <- data.table(
        ID=ID, AFT=AFT, AFN=AFN, ADT=ADT, DPT=DPT, ADN=ADN, DPN=DPN, ADT_FWD=adt_fwd, ADT_REV=adt_rev,
        str=str, str_ru=str_ru, str_len=str_len, str_diff=str_diff, ecnt=ecnt,
        tlod=tlod, nlod=nlod, flod=flod, xfet=xfet, sor=sor, mqrs=mqrs, mqs=mqs,
        oxog=oxog, pid=pid, c3_adx=c3_adx
    )
    return(var.qnt)
}

.pairBND <- function(bnd) {
    tmp1 <- data.table(chr1=as.character(seqnames(bnd)), pos1=start(bnd), ref1=as.character(ref(bnd)))
    tmp2 <- str_match(rowRanges(bnd)$ALT, "(.*)(\\[|])(.*):(.*)(\\[|])(.*)")
    tmp2 <- as.data.table(tmp2)
    tmp3 <- cbind(tmp1, tmp2)
    tmp3[V2!="", bnd1:="]"]
    tmp3[V7!="", bnd1:="["]
    tmp3[bnd1=="]", insert:=str_sub(V2, 2, -1)]
    tmp3[bnd1=="[", insert:=str_sub(V7, 1, -2)]
    tmp4 <- tmp3[,.(chr1, pos1, bnd1, chr2=V4, pos2=as.integer(V5), bnd2=V6, insert)]
    tmp5 <- data.table(
        id1=names(bnd),
        idx1=seq_len(length(bnd)),
        id2=info(bnd)$MATEID
    )
    setkey(tmp5, id1)
    tmp5$idx2 <- tmp5[J(tmp5$id2)]$idx1
    setkey(tmp5, idx1)
    bnd.pair <- cbind(tmp4, tmp5)
    bnd.pair <- bnd.pair[complete.cases(bnd.pair)] ## tnscope bug (maybe not dropped 'decoy' chroms?)
    bnd.pair[,bnd.id:=paste(
                  ifelse(idx1<idx2, idx1, idx2),
                  ifelse(idx1<idx2, idx2, idx1),
                  sep=":")]
    bnd.pair[,":="(
        chr1=factor(chr1, seqlevels(bnd), ordered=TRUE),
        chr2=factor(chr2, seqlevels(bnd), ordered=TRUE)
    )]
    bnd.pair <- bnd.pair[,.SD[1], bnd.id]
    bnd.pair[, proximal:=(chr1==chr2 & abs(pos1-pos2)<1e6)]
    setkey(bnd.pair, chr1, pos1, bnd1, chr2, pos2, bnd2)
    return(bnd.pair)
}

makeReport <- function(var, transcripts) {
    if ("SVTYPE" %in% names(info(var))) {
        var <- var[is.na(info(var)$SVTYPE)]
    }
    if (length(var)>0) {
        var.rng <- .varRNG(var)
        if (!all(unlist(is.na(info(var)$CSQ)))) {
            var.csq <- .varCSQ(var, transcripts)[,-1]
        } else {
            var.csq <- data.table(
                SYMBOL = NA_character_, TRANSCRIPT = NA_character_,
                EXON = NA_character_, HGVSc = NA_character_, HGVSp = NA_character_,
                Consequence = NA_character_, IMPACT = NA_character_, SIFT = NA_character_,
                Warnings = NA_character_, HGVS_Equivalent = NA_character_
            )
            var.csq <- var.csq[1:length(var)]
        }
        var.qnt <- .varQNT(var)[,-1]
        var.ann <- .varANN(var)[,-1]
        rep <- data.table(var.rng, var.csq, var.qnt, var.ann)
        for (i in names(which(sapply(rep, class)=="character"))) {
            rep[get(i)=="",(i):=NA_character_]
        }
    } else {
        ## create empty table
        rep <- structure(list(ID = character(0), chr = character(0), pos = integer(0), ref = character(0),
                              alt = character(0), filter = character(0), 
                              context = character(0), sbs96 = character(0), SYMBOL = character(0), GENE = character(0),
                              TRANSCRIPT = character(0), EXON = character(0), HGVSc = character(0), HGVSp = character(0),
                              Consequence = character(0), IMPACT = character(0), SIFT = character(0),
                              Warnings = character(0), HGVS_Equivalent = character(0),
                              AFT = numeric(0), AFN = numeric(0), ADT = numeric(0), DPT = numeric(0),
                              ADN = numeric(0), DPN = numeric(0), ADT_FWD = integer(0), ADT_REV = integer(0),
                              str = logical(0), str_ru=character(0), str_len=numeric(0), str_diff=numeric(0),
                              ecnt = numeric(0), tlod = numeric(0), nlod = numeric(0), flod = numeric(0), xfet = numeric(0),
                              sor = numeric(0), mqrs = numeric(0), mqs=numeric(0), oxog = numeric(0), pid = numeric(0),
                              c3_adx = character(0), dbsnp = integer(0), clinid = integer(0), clinvar = character(0), cosmic_cnt = integer(0),
                              gnomad_af = numeric(0), kg_af = numeric(0), pon = list(), rmsk_hit = list()), row.names = c(NA, 0L),
                              class = c("data.table", "data.frame"))
    }
    return(rep)
}

makeStructuralReport <- function(var, transcripts) {
    if ("SVTYPE" %in% names(info(var))) {
        var <- var[!is.na(info(var)$SVTYPE)]
    }
    if (nrow(var) != 0) {
        ## pair BND
        bnd.pair <- .pairBND(var)
        ## consequences for breakends
        csq.pick <- .pickCSQ2(var, bnd.pair, transcripts)
        ## consequences for breakend pairs
        csq.pick.short.1 <- csq.pick[bnd.pair$id1,.(
                                                       SYMBOL1=SYMBOL, GENE1=Gene, TRANSCRIPT1=id.short, EXON1=EXON, INTRON1=INTRON,
                                                       STRAND1=factor(ifelse(STRAND=="1", "+", ifelse(STRAND=="-1", "-", "*")),
                                                       levels=c("+", "-", "*"), ordered=TRUE)
                                                   )]
        csq.pick.short.2 <- csq.pick[bnd.pair$id2,.(
                                                       SYMBOL2=SYMBOL, GENE2=Gene, TRANSCRIPT2=id.short, EXON2=EXON, INTRON2=INTRON,
                                                       STRAND2=factor(ifelse(STRAND=="1", "+", ifelse(STRAND=="-1", "-", "*")),
                                                       levels=c("+", "-", "*"), ordered=TRUE)
                                                   )]
        ##  coverage for breakends
        adt <- geno(var)$AD[bnd.pair$idx1,1,2] + geno(var)$AD[bnd.pair$idx2,1,2]
        dpt <- rowSums(geno(var)$AD[bnd.pair$idx1,1,,drop=FALSE] + geno(var)$AD[bnd.pair$idx2,1,,drop=FALSE])
        qual <- rowRanges(var)$QUAL[bnd.pair$idx1]

        ## merge breakends with consequences and coverage
        bnd.csq <- data.table(
            bnd.pair,
            AFT=adt/dpt,
            ADT=adt,
            DPT=dpt,
            QUAL=qual,
            MANTA.FILTER=info(var)$MANTA_FILTER[bnd.pair$idx1],
            MANTA.SCORE=as.integer(info(var)$MANTA_SOMATICSCORE[bnd.pair$idx1]),
            MANTA.CONTIG=info(var)$MANTA_CONTIG[bnd.pair$idx1],
            csq.pick.short.1,
            csq.pick.short.2
        )

        ## basic classification
        bnd.csq[is.na(STRAND1), STRAND1:="*"]
        bnd.csq[is.na(STRAND2), STRAND2:="*"]
        bnd.csq[,topo:="inversion"]
        bnd.csq[chr1!=chr2, topo:="translocation"]
        bnd.csq[chr1==chr2 & bnd1=="[" & bnd2=="]", topo:="duplication"]
        bnd.csq[chr1==chr2 & bnd1=="]" & bnd2=="[", topo:="deletion"]
        bnd.csq[(bnd1=="[" & STRAND1=="+") | (bnd1=="]" & STRAND1=="-"), END1:="3p"]
        bnd.csq[(bnd1=="]" & STRAND1=="+") | (bnd1=="[" & STRAND1=="-"), END1:="5p"]
        bnd.csq[(bnd2=="[" & STRAND2=="+") | (bnd2=="]" & STRAND2=="-"), END2:="3p"]
        bnd.csq[(bnd2=="]" & STRAND2=="+") | (bnd2=="[" & STRAND2=="-"), END2:="5p"]

        ## flip strands
        flip <- bnd.csq[,
            ifelse(
                SYMBOL1==SYMBOL2 & STRAND1=="-", TRUE,
            ifelse(
                SYMBOL1==SYMBOL2 & STRAND1=="+", FALSE,
            ifelse(
                !is.na(END1) & !is.na(END2) & END1=="3p" & END2=="5p", TRUE,
            ifelse(
                !is.na(END1) & is.na(END2) & END1=="3p", TRUE,
            ifelse(
                !is.na(END2) & is.na(END1) & END2=="5p", TRUE,
                FALSE
            )))))]
        bnd.csq.p <- bnd.csq[!flip,
                             .(bnd.id, chr1, pos1, chr2, pos2, topo, insert, QUAL, AFT, ADT, DPT,
                               MANTA.FILTER, MANTA.SCORE, MANTA.CONTIG,
                               SYMBOL1,GENE1,TRANSCRIPT1,EXON1,INTRON1,STRAND1,END1,
                               SYMBOL2,GENE2,TRANSCRIPT2,EXON2,INTRON2,STRAND2,END2)]
        bnd.csq.m <- bnd.csq[flip,
                             .(bnd.id, chr2, pos2, chr1, pos1, topo, insert, QUAL, AFT, ADT, DPT,
                               MANTA.FILTER, MANTA.SCORE, MANTA.CONTIG,
                               SYMBOL2,GENE2,TRANSCRIPT2,EXON2,INTRON2,STRAND2,END2,
                               SYMBOL1,GENE1,TRANSCRIPT1,EXON1,INTRON1,STRAND1,END1)]
        rep <- rbind(bnd.csq.p, unname(bnd.csq.m))
        for (i in names(which(sapply(rep, class)=="character"))) {
            rep[get(i)=="",(i):=NA_character_]
        }
    } else {
        ## create empty table
        rep <- data.table(1)[,`:=`(c(
                             "bnd.id", "chr1", "pos1", "chr2", "pos2", "topo", "insert",
                             "QUAL", "AFT", "ADT", "DPT", "MANTA.FILTER", "MANTA.SCORE", "MANTA.CONTIG",
                             "SYMBOL1", "GENE1", "TRANSCRIPT1", "EXON1", "INTRON1", "STRAND1", "END1",
                             "SYMBOL2", "GENE2", "TRANSCRIPT2", "EXON2", "INTRON2", "STRAND2", "END2"),NA)
                             ][,V1:=NULL][.0]
    }
    return(rep)
}
