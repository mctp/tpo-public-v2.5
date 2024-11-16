#!/usr/bin/env bash

function merge-crisp-align {
    local ALN=$1
    local INPUT=$2
    local PFX=$3
    local ALIGN_FASTA=$4
    local ALIG_CRAM=''
    local ALIG_SJ=''
    local CHIM_PE_CRAM=''
    local CHIM_SE_CRAM=''
    local CHIM_PE_JCT=''
    local CHIM_SE_JCT=''
    local FQ1=''
    local FQ2=''
    IFS=';' read -ra ALNA <<< $ALN
    for i in "${ALNA[@]}"; do
        local ID=$(basename $i)
        ALIG_CRAM+=" $INPUT/$ID/$ID-alig.cram"
        CHIM_PE_CRAM+=" $INPUT/$ID/$ID-chim-pe.cram"
        CHIM_SE_CRAM+=" $INPUT/$ID/$ID-chim-se.cram"
        ALIG_SJ+=" $INPUT/$ID/$ID-sj.tab.gz"
        CHIM_PE_JCT+=" $INPUT/$ID/$ID-chim-pe.jnc.gz"
        CHIM_SE_JCT+=" $INPUT/$ID/$ID-chim-se.jnc.gz"
        FQ1+=" $INPUT/$ID/$ID-unmapped_1.fq.gz"
        FQ2+=" $INPUT/$ID/$ID-unmapped_2.fq.gz"
    done
    local ALNN="${#ALNA[@]}"

    samtools merge -f -@ 8 -O BAM --reference $ALIGN_FASTA $PFX-alig.bam $ALIG_CRAM &
    samtools merge -f -@ 4 -O BAM --reference $ALIGN_FASTA $PFX-chim-pe.bam $CHIM_PE_CRAM &
    samtools merge -f -@ 4 -O BAM --reference $ALIGN_FASTA $PFX-chim-se.bam $CHIM_SE_CRAM &
    wait
    samtools index -@ 8 $PFX-alig.bam
    samtools index -@ 8 $PFX-chim-pe.bam
    samtools index -@ 8 $PFX-chim-se.bam
    
    if (( ALNN > 1 )); then
        cat $ALIG_SJ > $PFX-sj.tab.gz
        cat $CHIM_PE_JCT > $PFX-chim-pe.jnc.gz
        cat $CHIM_SE_JCT > $PFX-chim-se.jnc.gz
        cat $FQ1 > $PFX-unmapped_1.fq.gz
        cat $FQ2 > $PFX-unmapped_2.fq.gz
    elif [ "$INPUT/$ID/$ID" != "$PFX" ]; then
        ln $ALIG_SJ $PFX-sj.tab.gz
        ln $CHIM_PE_JCT $PFX-chim-pe.jnc.gz
        ln $CHIM_SE_JCT $PFX-chim-se.jnc.gz
        ln $FQ1 $PFX-unmapped_1.fq.gz
        ln $FQ2 $PFX-unmapped_2.fq.gz
    fi

}
