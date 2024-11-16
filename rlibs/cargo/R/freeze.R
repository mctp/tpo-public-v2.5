#' @export
fusionFreeze <- function(mvf) {
    tmp <- svFormat(mvf$tables$fusion, FALSE)
    cbind(id=mvf$tables$fusion$id, tmp)
}



expression2df <- function(TBL, VALUE="rpkm", ROUND=NA){
  if(!is.na(ROUND)){
    TBL[[VALUE]] <- round(TBL[[VALUE]], ROUND)
  }
  TBL.df <- reshape2::dcast(TBL, GENE~id, value.var=VALUE)
  TBL.df <- tibble::add_column(TBL.df, SYMBOL=TBL$SYMBOL[match(TBL.df$GENE,TBL$GENE)], .after=1)
  return(TBL.df)
}




somatic2Maf <- function(TBL) {
  TBL[["ADT"]] <- as.numeric(TBL[["ADT"]])
  TBL[["ADN"]] <- as.numeric(TBL[["ADN"]])
  TBL[["AFT"]] <- as.numeric(TBL[["AFT"]])
  TBL[["AFN"]] <- as.numeric(TBL[["AFN"]])
    maf <- TBL[,.(
        Patient_ID = id, #not in standard MAF, but should be included?
        Hugo_Symbol=gene_name,
        Entrez_Gene_Id=".",
        Center="MCTP",
        NCBI_Build="GRCh38",
        Chromosome=chr,
        Start_Position=
            ifelse(str_length(ref)<=str_length(alt), pos, # SNV, MNV, insertion
            ifelse(str_length(ref)>str_length(alt) & str_sub(ref,1,1)==str_sub(alt,1,1), pos + 1, # deletion
            ifelse(str_length(ref)>str_length(alt) & str_sub(ref,1,1)!=str_sub(alt,1,1), pos, # deletion
                   -1)))
       ,
        End_Position=pos + str_length(ref) - 1,
        Strand="+",
        Variant_Classification=
            ifelse(
                grepl("start_lost", Consequence), "Translation_Start_Site",
            ifelse(
                grepl("stop_gained", Consequence), "Nonsense_Mutation",
            ifelse(
                grepl("stop_lost", Consequence), "Nonstop_Mutation",
            ifelse(
                grepl("splice_donor_variant|splice_acceptor_variant", Consequence), "Splice_Site",
            ifelse(
                grepl("frameshift_variant", Consequence) & str_length(ref)<str_length(alt), "Frame_Shift_Ins",
            ifelse(
                grepl("frameshift_variant", Consequence) & str_length(ref)>str_length(alt), "Frame_Shift_Del",
            ifelse(
                grepl("inframe_insertion", Consequence), "In_Frame_Ins",
            ifelse(
                grepl("inframe_deletion", Consequence), "In_Frame_Del",
            ifelse(
                grepl("missense_variant", Consequence), "Missense_Mutation",
            ifelse(
                grepl("protein_altering_variant", Consequence), "Protein_Altering",
            ifelse(
                grepl("upstream_gene_variant", Consequence), "Missense_Mutation",
            "Silent"))))))))))),
        Variant_Type=
            ifelse(str_length(ref)==str_length(alt) & str_length(ref)==1, "SNP",
            ifelse(str_length(ref)==str_length(alt) & str_length(ref)==2, "DNP",
            ifelse(str_length(ref)==str_length(alt) & str_length(ref)==3, "TNP",
            ifelse(str_length(ref)==str_length(alt) & str_length(ref)>=4, "ONP",
            ifelse(str_length(ref)<str_length(alt), "INS",
            ifelse(str_length(ref)>str_length(alt), "DEL",
                   "XXX"
                   ))))))
       ,
        Reference_Allele=
            ifelse(str_length(ref)==str_length(alt), ref,
            ifelse(str_length(ref)<str_length(alt) & str_length(ref)==1, "-",
            ifelse(str_length(ref)>str_length(alt) & str_length(alt)==1, str_sub(ref, 2),
            ifelse(str_sub(ref,1,1)==str_sub(alt,1,1), str_sub(ref, 2),
                   ref
                   )))),
        Tumor_Seq_Allele1=
            ifelse(str_length(ref)==str_length(alt), ref,
            ifelse(str_length(ref)<str_length(alt) & str_length(ref)==1, "-",
            ifelse(str_length(ref)>str_length(alt) & str_length(alt)==1, str_sub(ref, 2),
            ifelse(str_sub(ref,1,1)==str_sub(alt,1,1), str_sub(ref, 2),
                   ref
                   )))),
        Tumor_Seq_Allele2=
            ifelse(str_length(ref)==str_length(alt), alt,
            ifelse(str_length(ref)<str_length(alt) & str_length(ref)==1, str_sub(alt, 2), 
            ifelse(str_length(ref)>str_length(alt) & str_length(alt)==1, "-",
            ifelse(str_sub(ref,1,1)==str_sub(alt,1,1), str_sub(alt, 2),
                   alt
                   )))),
        dbSNP_RS=dbsnp,
        dbSNP_Val_Status="by1000genomes",
        Tumor_Sample_Barcode=sapply(str_split(id, "\\."), "[", 2),
        Matched_Norm_Sample_Barcode=sapply(str_split(id, "\\."), "[", 3),
        Match_Norm_Seq_Allele1=
            ifelse(str_length(ref)==str_length(alt), ref,
            ifelse(str_length(ref)<str_length(alt) & str_length(ref)==1, "-",
            ifelse(str_length(ref)>str_length(alt) & str_length(alt)==1, str_sub(ref, 2),
            ifelse(str_sub(ref,1,1)==str_sub(alt,1,1), str_sub(ref, 2),
                   ref
                   )))),
        Match_Norm_Seq_Allele2=
            ifelse(str_length(ref)==str_length(alt), ref,
            ifelse(str_length(ref)<str_length(alt) & str_length(ref)==1, "-",
            ifelse(str_length(ref)>str_length(alt) & str_length(alt)==1, str_sub(ref, 2),
            ifelse(str_sub(ref,1,1)==str_sub(alt,1,1), str_sub(ref, 2),
                   ref
                   )))),
        Tumor_Validation_Allele1=".",
        Tumor_Validation_Allele2=".",
        Match_Norm_Validation_Allele1=".",
        Match_Norm_Validation_Allele2=".",
        Verification_Status="Unknown",
        Validation_Status="Untested",
        Mutation_Status="Somatic",
        Sequencing_Phase=".",
        Sequence_Source="WXS",
        Validation_Method="none",
        Score=".",
        BAM_File=".",
        Sequencer="Illumina HiSeq",
        Tumor_Sample_UUID=".",
        Normal_Sample_UUID=".",
        HGVSc=HGVSc,
        HGVSp=HGVSp,
        HGVSp_Short=Add_HGVSp_short(HGVSp),
        Transcript_ID=TRANSCRIPT,
        Exon_Number=EXON,
        t_depth=DPT,
        t_ref_count=DPT-ADT,
        t_alt_count=ADT,
        n_depth=DPN,
        n_ref_count=DPN-ADN,
        n_alt_count=ADN,
       # all_effects
       # Allele
       # Gene
       # Feature
       # Feature_type
       # One_Consequence
       # Consequence
       # cDNA_position
       # CDS_position
       # Protein_position
       # Amino_acids
       # Codons
       # Existing_variation
       # ALLELE_NUM
       # DISTANCE #this is #60 in the list
        IMPACT=IMPACT, #94
        SIFT=SIFT, #73
        GMAF=kg_af,
        COSMIC=cosmic_cnt,
        CLIN_SIG=clinvar,
        ExAC_AF=gnomad_af,
        pon=pon
    )
    ]
    return(maf)
}

## 1 - Hugo_Symbol	  HUGO symbol for the gene (HUGO symbols are always in all caps). "Unknown" is used for regions that do not correspond to a gene
## 2 - Entrez_Gene_Id	  Entrez gene ID (an integer). "0" is used for regions that do not correspond to a gene region or Ensembl ID
## 3 - Center	One or more genome sequencing center reporting the variant
## 4 - NCBI_Build	The reference genome used for the alignment (GRCh38)
## 5 - Chromosome	The affected chromosome (chr1)
## 6 - Start_Position	Lowest numeric position of the reported variant on the genomic reference sequence. Mutation start coordinate
## 7 - End_Position	Highest numeric genomic position of the reported variant on the genomic reference sequence. Mutation end coordinate
## 8 - Strand	Genomic strand of the reported allele. Currently, all variants will report the positive strand: '+'
## 9 - Variant_Classification	Translational effect of variant allele
## 10 - Variant_Type	Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP (di-nucleotide polymorphism) but for three consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of four or more (SNP, DNP, TNP, ONP, INS, DEL, or Consolidated)
## 11 - Reference_Allele	The plus strand reference allele at this position. Includes the deleted sequence for a deletion or "-" for an insertion
## 12 - Tumor_Seq_Allele1	Primary data genotype for tumor sequencing (discovery) allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases
## 13 - Tumor_Seq_Allele2	Tumor sequencing (discovery) allele 2
## 14 - dbSNP_RS	The rs-IDs from the   dbSNP database, "novel" if not found in any database used, or null if there is no dbSNP record, but it is found in other databases
## 15 - dbSNP_Val_Status	The dbSNP validation status is reported as a semicolon-separated list of statuses. The union of all rs-IDs is taken when there are multiple
## 16 - Tumor_Sample_Barcode	Aliquot barcode for the tumor sample
## 17 - Matched_Norm_Sample_Barcode	Aliquot barcode for the matched normal sample
## 18 - Match_Norm_Seq_Allele1	Primary data genotype. Matched normal sequencing allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases (cleared in somatic MAF)
## 19 - Match_Norm_Seq_Allele2	Matched normal sequencing allele 2
## 20 - Tumor_Validation_Allele1	Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases
## 21 - Tumor_Validation_Allele2	Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 2
## 22 - Match_Norm_Validation_Allele1	Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases (cleared in somatic MAF)
## 23 - Match_Norm_Validation_Allele2	Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 2 (cleared in somatic MAF)
## 24 - Verification_Status	Second pass results from independent attempt using same methods as primary data source. Generally reserved for 3730 Sanger Sequencing
## 25 - Validation_Status	Second pass results from orthogonal technology
## 26 - Mutation_Status	An assessment of the mutation as somatic, germline, LOH, post transcriptional modification, unknown, or none. The values allowed in this field are constrained by the value in the Validation_Status field
## 27 - Sequencing_Phase	TCGA sequencing phase (if applicable). Phase should change under any circumstance that the targets under consideration change
## 28 - Sequence_Source	Molecular assay type used to produce the analytes used for sequencing. Allowed values are a subset of the SRA 1.5 library_strategy field values. This subset matches those used at CGHub
## 29 - Validation_Method	The assay platforms used for the validation call
## 30 - Score	Not in use
## 31 - BAM_File	Not in use
## 32 - Sequencer	Instrument used to produce primary sequence data
## 33 - Tumor_Sample_UUID	GDC aliquot UUID for tumor sample
## 34 - Matched_Norm_Sample_UUID	GDC aliquot UUID for matched normal sample
## 35 - HGVSc	The coding sequence of the variant in HGVS recommended format
## 36 - HGVSp	The protein sequence of the variant in HGVS recommended format. "p.=" signifies no change in the protein
## 37 - HGVSp_Short	Same as the HGVSp column, but using 1-letter amino-acid codes
## 38 - Transcript_ID	  Ensembl ID of the transcript affected by the variant
## 39 - Exon_Number	The exon number (out of total number)
## 40 - t_depth	Read depth across this locus in tumor BAM
## 41 - t_ref_count	Read depth supporting the reference allele in tumor BAM
## 42 - t_alt_count	Read depth supporting the variant allele in tumor BAM
## 43 - n_depth	Read depth across this locus in normal BAM
## 44 - n_ref_count	Read depth supporting the reference allele in normal BAM (cleared in somatic MAF)
## 45 - n_alt_count	Read depth supporting the variant allele in normal BAM (cleared in somatic MAF)
## 46 - all_effects	A semicolon delimited list of all possible variant effects, sorted by priority ([Symbol,Consequence,HGVSp_Short,Transcript_ID,RefSeq,HGVSc,Impact,Canonical,Sift,PolyPhen,Strand])
## 47 - Allele	The variant allele used to calculate the consequence
## 48 - Gene	Stable Ensembl ID of affected gene
## 49 - Feature	Stable Ensembl ID of feature (transcript, regulatory, motif)
## 50 - Feature_type	Type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature (or blank)
## 51 - One_Consequence	The single consequence of the canonical transcript in   sequence ontology terms
## 52 - Consequence	Consequence type of this variant;   sequence ontology terms
## 53 - cDNA_position	Relative position of base pair in the cDNA sequence as a fraction. A "-" symbol is displayed as the numerator if the variant does not appear in cDNA
## 54 - CDS_position	Relative position of base pair in coding sequence. A "-" symbol is displayed as the numerator if the variant does not appear in coding sequence
## 55 - Protein_position	Relative position of affected amino acid in protein. A "-" symbol is displayed as the numerator if the variant does not appear in coding sequence
## 56 - Amino_acids	Only given if the variation affects the protein-coding sequence
## 57 - Codons	The alternative codons with the variant base in upper case
## 58 - Existing_variation	Known identifier of existing variation
## 59 - ALLELE_NUM	Allele number from input; 0 is reference, 1 is first alternate etc.
## 60 - DISTANCE	Shortest distance from the variant to transcript
## 61 - TRANSCRIPT_STRAND	The DNA strand (1 or -1) on which the transcript/feature lies
## 62 - SYMBOL	The gene symbol
## 63 - SYMBOL_SOURCE	The source of the gene symbol
## 64 - HGNC_ID	Gene identifier from the HUGO Gene Nomenclature Committee if applicable
## 65 - BIOTYPE	Biotype of transcript
## 66 - CANONICAL	A flag (YES) indicating that the VEP-based canonical transcript, the longest translation, was used for this gene. If not, the value is null
## 67 - CCDS	The   CCDS identifier for this transcript, where applicable
## 68 - ENSP	The Ensembl protein identifier of the affected transcript
## 69 - SWISSPROT	  UniProtKB/Swiss-Prot accession
## 70 - TREMBL	UniProtKB/TrEMBL identifier of protein product
## 71 - UNIPARC	UniParc identifier of protein product
## 72 - RefSeq	RefSeq identifier for this transcript
## 73 - SIFT	The   SIFT prediction and/or score, with both given as prediction (score)
## 74 - PolyPhen	The   PolyPhen prediction and/or score
## 75 - EXON	The exon number (out of total number)
## 76 - INTRON	The intron number (out of total number)
## 77 - DOMAINS	The source and identifier of any overlapping protein domains
## 78 - GMAF	Non-reference allele and frequency of existing variant in   1000 Genomes
## 79 - AFR_MAF	Non-reference allele and frequency of existing variant in 1000 Genomes combined African population
## 80 - AMR_MAF	Non-reference allele and frequency of existing variant in 1000 Genomes combined American population
## 81 - ASN_MAF	Non-reference allele and frequency of existing variant in 1000 Genomes combined Asian population
## 82 - EAS_MAF	Non-reference allele and frequency of existing variant in 1000 Genomes combined East Asian population
## 83 - EUR_MAF	Non-reference allele and frequency of existing variant in 1000 Genomes combined European population
## 84 - SAS_MAF	Non-reference allele and frequency of existing variant in 1000 Genomes combined South Asian population
## 85 - AA_MAF	Non-reference allele and frequency of existing variant in   NHLBI-ESP African American population
## 86 - EA_MAF	Non-reference allele and frequency of existing variant in NHLBI-ESP European American population
## 87 - CLIN_SIG	Clinical significance of variant from dbSNP
## 88 - SOMATIC	Somatic status of each ID reported under Existing_variation (0, 1, or null)
## 89 - PUBMED	Pubmed ID(s) of publications that cite existing variant
## 90 - MOTIF_NAME	The source and identifier of a transcription factor binding profile aligned at this position
## 91 - MOTIF_POS	The relative position of the variation in the aligned TFBP
## 92 - HIGH_INF_POS	A flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP) (Y, N, or null)
## 93 - MOTIF_SCORE_CHANGE	The difference in motif score of the reference and variant sequences for the TFBP
## 94 - IMPACT	The impact modifier for the consequence type
## 95 - PICK	Indicates if this block of consequence data was picked by VEP's   pick feature (1 or null)
## 96 - VARIANT_CLASS	Sequence Ontology variant class
## 97 - TSL	  Transcript support level, which is based on independent RNA analyses
## 98 - HGVS_OFFSET	Indicates by how many bases the HGVS notations for this variant have been shifted
## 99 - PHENO	Indicates if existing variant is associated with a phenotype, disease or trait (0, 1, or null)
## 100 - MINIMISED	Alleles in this variant have been converted to minimal representation before consequence calculation (1 or null)
## 101 - ExAC_AF	Global Allele Frequency from   ExAC
## 102 - ExAC_AF_Adj	Adjusted Global Allele Frequency from ExAC
## 103 - ExAC_AF_AFR	African/African American Allele Frequency from ExAC
## 104 - ExAC_AF_AMR	American Allele Frequency from ExAC
## 105 - ExAC_AF_EAS	East Asian Allele Frequency from ExAC
## 106 - ExAC_AF_FIN	Finnish Allele Frequency from ExAC
## 107 - ExAC_AF_NFE	Non-Finnish European Allele Frequency from ExAC
## 108 - ExAC_AF_OTH	Other Allele Frequency from ExAC
## 109 - ExAC_AF_SAS	South Asian Allele Frequency from ExAC
## 110 - GENE_PHENO	Indicates if gene that the variant maps to is associated with a phenotype, disease or trait (0, 1, or null)
## 111 - FILTER	Copied from input VCF. This includes filters implemented directly by the variant caller and other external software used in the DNA-Seq pipeline. See below for additional details.
## 112 - CONTEXT	The reference allele per VCF specs, and its five flanking base pairs
## 113 - src_vcf_id	GDC UUID for the input VCF file
## 114 - tumor_bam_uuid	GDC UUID for the tumor bam file
## 115 - normal_bam_uuid	GDC UUID for the normal bam file
## 116 - case_id	GDC UUID for the case
## 117 - GDC_FILTER	GDC filters applied universally across all MAFs
## 118 - COSMIC	Overlapping COSMIC variants
## 119 - MC3_Overlap	Indicates whether this region overlaps with an MC3 variant for the same sample pair
## 120 - GDC_Validation_Status	GDC implementation of validation checks. See notes section (#5) below for details
## 121 - GDC_Valid_Somatic	True or False (not in somatic MAF)
## 122 - vcf_region	Colon separated string containing the CHROM, POS, ID, REF, and ALT columns from the VCF file (e.g., chrZ:20:rs1234:A:T) (not in somatic MAF)
## 123 - vcf_info	INFO column from VCF (not in somatic MAF)
## 124 - vcf_format	FORMAT column from VCF (not in somatic MAF)
## 125 - vcf_tumor_gt	Tumor sample genotype column from VCF (not in somatic MAF)
## 126 - vcf_normal_gt	Normal sample genotype column from VCF (not in somatic MAF)




Add_HGVSp_short <- function(HGVSp){
  proteinlib <- c('Cys'= 'C', 'Asp'= 'D', 'Ser'= 'S', 'Gln'= 'Q', 'Lys'= 'K',
                  'Ile'= 'I', 'Pro'= 'P', 'Thr'= 'T', 'Phe'= 'F', 'Asn'= 'N', 
                  'Gly'= 'G', 'His'= 'H', 'Leu'= 'L', 'Arg'= 'R', 'Trp'= 'W', 
                  'Ala'= 'A', 'Val'='V', 'Glu'= 'E', 'Tyr'= 'Y', 'Met'= 'M')
  HGVSp <- gsub("p[.]", "", HGVSp)
  for(i in 1:length(HGVSp)){
    for(p in names(proteinlib)){
      if( grepl(p, HGVSp[i], ignore.case = T) ){
        HGVSp[i] <- gsub(p, proteinlib[p], HGVSp[i])
      }
    }
  }
  HGVSp[is.na(HGVSp)] <- ""
  return(HGVSp)
}












