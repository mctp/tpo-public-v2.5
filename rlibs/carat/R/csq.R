#' @export
structuralCSQ <- function(tbl) {
    ## compute consequences
    tbl$CSQ1 <- NA_character_
    tbl$CSQ2 <- NA_character_
    tbl[!is.na(GENE1) & END1=="5p", ":="(CSQ1="truncation")]
    tbl[!is.na(GENE1) & END1=="3p", ":="(CSQ1="start-loss")]
    tbl[!is.na(GENE2) & END2=="5p", ":="(CSQ2="truncation")]
    tbl[!is.na(GENE2) & END2=="3p", ":="(CSQ2="start-loss")]
    tbl[GENE1==GENE2 & topo=="duplication", ":="(CSQ1="internal-duplication", CSQ2="internal-duplication")]
    tbl[GENE1==GENE2 & topo=="deletion", ":="(CSQ1="internal-deletion", CSQ2="internal-deletion")]
    tbl[GENE1==GENE2 & topo=="inversion", ":="(CSQ1="internal-inversion", CSQ2="internal-inversion")]
    tbl[GENE1!=GENE2 & END1=="5p" & END2=="3p", ":="(CSQ1="fusion5", CSQ2="fusion3")]
    ## compute affected exons
    ex1 <- as.integer(str_match(tbl$EXON1, "([0-9]*)/")[,2])
    in1 <- as.integer(str_match(tbl$INTRON1, "([0-9]*)/")[,2])
    ex2 <- as.integer(str_match(tbl$EXON2, "([0-9]*)/")[,2])
    in2 <- as.integer(str_match(tbl$INTRON2, "([0-9]*)/")[,2])
    b1 <- ifelse(!is.na(ex1), ex1, ifelse(tbl$END1=="5p", in1+1, in1))
    b2 <- ifelse(!is.na(ex2), ex2, ifelse(tbl$END2=="5p", in2+1, in2))
    ## crazy logic for internal
    b1[is.na(b1)] <- ""
    b2 <- ifelse(!is.na(ex2), ex2, in2)
    b2[is.na(b2)] <- ""
    tbl$B1 <- b1
    tbl$B2 <- b2
    tbl$B12 <- ifelse(b1==b2, b1,
               ifelse(b1=="", paste0("--", b2),
               ifelse(b2=="", paste0(b1, "--"),
               ifelse(as.integer(b1)<=as.integer(b2), paste(b1, b2, sep="--"), NA))))
    ## set exons
    tbl$EXONS1 <- NA_character_
    tbl$EXONS2 <- NA_character_
    tbl[CSQ1 %in% c("truncation", "fusion5") & !is.na(B1), EXONS1:=paste0(B1, "--")]
    tbl[CSQ1 %in% c("start-loss", "fusion3") & !is.na(B1), EXONS1:=paste0("--", B1)]
    tbl[CSQ2 %in% c("truncation", "fusion5") & !is.na(B2), EXONS2:=paste0(B2, "--")]
    tbl[CSQ2 %in% c("start-loss", "fusion3") & !is.na(B2), EXONS2:=paste0("--", B2)]
    tbl[CSQ1 %in% c("internal-deletion", "internal-duplication"), ":="(EXONS1=B12, EXONS2=B12)]
    structural.tbl <- tbl
    return(structural.tbl)
}

#' @export
fusionCSQ <- function(tbl) {
    tbl$CSQF5 <- NA_character_
    tbl$CSQF3 <- NA_character_
    tbl$CSQR5 <- NA_character_
    tbl$CSQR3 <- NA_character_
    tbl[(gene_names.5.1 != "" & gene_names.3.1 != ""),
        ":="(CSQF5="fusion5", CSQF3="fusion3")
        ]
    tbl[(gene_names.5.1 != "" & gene_names.3.1 == ""),
        ":="(CSQF5="truncation")
        ]
    tbl[(gene_names.5.1 == "" & gene_names.3.1 != ""),
        ":="(CSQF3="start-loss")
        ]
    tbl[(gene_names.5.2 != ""),
        ":="(CSQR5="truncation")
        ]
    tbl[(gene_names.3.2 != ""),
        ":="(CSQR3="truncation")
        ]
    fusion.tbl <- tbl
    return(fusion.tbl)
}
