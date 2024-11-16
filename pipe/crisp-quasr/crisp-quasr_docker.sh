#!/usr/bin/env bash
echo $(date) >> /job/docker_start

set -eu
shopt -s nullglob

#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate crisp

source /code/pipe/crisp-align/crisp-align_helper.sh

source /job/config.txt

NCORES=$(nproc --all)

mkdir -p /input/$ID

PFX=/output/$ID

merge-crisp-align $ALN /input /input/$ID/$ID $ALIGN_FASTA

if [ "$QUASR_KEEP" = true ]; then

    ln /job/input/$ID/$ID-alig.bam /job/output/$ID-alig.bam
    ln /job/input/$ID/$ID-alig.bam.bai /job/output/$ID-alig.bam.bai
    ln /job/input/$ID/$ID-chim-pe.bam /job/output/$ID-chim-pe.bam
    ln /job/input/$ID/$ID-chim-pe.bam.bai /job/output/$ID-chim-pe.bam.bai
    ln /job/input/$ID/$ID-chim-se.bam /job/output/$ID-chim-se.bam
    ln /job/input/$ID/$ID-chim-se.bam.bai /job/output/$ID-chim-se.bam.bai
    ln /job/input/$ID/$ID-sj.tab.gz /job/output/$ID-sj.tab.gz
    ln /job/input/$ID/$ID-chim-pe.jnc.gz /job/output/$ID-chim-pe.jnc.gz
    ln /job/input/$ID/$ID-chim-se.jnc.gz /job/output/$ID-chim-se.jnc.gz
    ln /job/input/$ID/$ID-unmapped_1.fq.gz /job/output/$ID-unmapped_1.fq.gz
    ln /job/input/$ID/$ID-unmapped_2.fq.gz /job/output/$ID-unmapped_2.fq.gz

fi

if [ "$QUASR_COUNT" = true ]; then

    featureCounts \
        --tmpDir /tmp \
        -T $NCORES \
        -a $QUASR_COUNTGTF \
        -o $PFX-count \
        $QUASR_COUNTARGS \
        /input/$ID/$ID-alig.bam &> $PFX-count.log
    
fi

if [ "$QUASR_MIXCR" = true ]; then
    
    source /code/pipe/crisp-quasr/$QUASR_MIXCRARGS.txt

    ## temparary files
    TRA=$(tempfile -d /tmp)
    TRB=$(tempfile -d /tmp)
    TRG=$(tempfile -d /tmp)
    TRX=$(tempfile -d /tmp)
    TRXN=$(tempfile -d /tmp)
    TRX_1=$(tempfile -d /tmp)
    TRX_2=$(tempfile -d /tmp)

    IGH=$(tempfile -d /tmp)
    IGL=$(tempfile -d /tmp)
    IGK=$(tempfile -d /tmp)
    IGX=$(tempfile -d /tmp)
    IGXN=$(tempfile -d /tmp)
    IGX_1=$(tempfile -d /tmp)
    IGX_2=$(tempfile -d /tmp)
    SEQ_1=$(tempfile -d /tmp -s fq.gz)
    SEQ_2=$(tempfile -d /tmp -s fq.gz)

    ## extract mapped reads
    samtools view -b /input/$ID/$ID-alig.bam $MIXCR_TRA > $TRA
    samtools view -b /input/$ID/$ID-alig.bam $MIXCR_TRB > $TRB
    samtools view -b /input/$ID/$ID-alig.bam $MIXCR_TRG > $TRG
    samtools merge -f -@ $NCORES $TRX $TRA $TRB $TRG
    samtools sort -n -@ $NCORES $TRX > $TRXN
    bedtools bamtofastq -i $TRXN -fq $TRX_1 -fq2 $TRX_2

    samtools view -b /input/$ID/$ID-alig.bam $MIXCR_IGH > $IGH
    samtools view -b /input/$ID/$ID-alig.bam $MIXCR_IGL > $IGL
    samtools view -b /input/$ID/$ID-alig.bam $MIXCR_IGK > $IGK
    samtools merge -f -@ $NCORES $IGX $IGH $IGL $IGK
    samtools sort -n -@ $NCORES $IGX > $IGXN
    bedtools bamtofastq -i $IGXN -fq $IGX_1 -fq2 $IGX_2

    ## compress
    pigz -p $NCORES $TRX_1
    pigz -p $NCORES $TRX_2
    pigz -p $NCORES $IGX_1
    pigz -p $NCORES $IGX_2

    cat $TRX_1.gz $IGX_1.gz /input/$ID/$ID-unmapped_1.fq.gz > $SEQ_1
    cat $TRX_2.gz $IGX_2.gz /input/$ID/$ID-unmapped_2.fq.gz > $SEQ_2

    ## cleanup
    rm $TRA $TRB $TRG $TRX $TRXN $IGH $IGL $IGK $IGX $IGXN $TRX_1.gz $TRX_2.gz $IGX_1.gz $IGX_2.gz

    ## MiXCR
    TEMP_3=$(tempfile -d /tmp)
    TEMP_3B=$(tempfile -d /tmp)
    TEMP_3T=$(tempfile -d /tmp)
    TEMP_3F=$(tempfile -d /tmp)
    TEMP_4=$(tempfile -d /tmp)
    TEMP_5=$(tempfile -d /tmp)

    conda deactivate ## MiXCR crashes with Conda OpenJDK (java.lang.ClassFormatError: Incompatible magic value 3893893500 in class file org/xml/sax/InputSource)
    mixcr align -f -t $NCORES -p rna-seq -s $MIXCR_SPECIES -OsaveOriginalReads=true -OallowPartialAlignments=true -r $PFX-mixcr-aln.rep.txt $SEQ_1 $SEQ_2 $TEMP_3
    mixcr filterAlignments -f -n 250000 -c IG  $TEMP_3 $TEMP_3B
    mixcr filterAlignments -f -n 250000 -c TCR $TEMP_3 $TEMP_3T
    rm $TEMP_3F
    mixcr mergeAlignments $TEMP_3B $TEMP_3T $TEMP_3F
    mixcr assemblePartial -f -r $PFX-mixcr-fix1.rep.txt $TEMP_3F $TEMP_4
    mixcr assemblePartial -f -r $PFX-mixcr-fix2.rep.txt $TEMP_4  $TEMP_5
    mixcr extend -f -r $PFX-mixcr-ext.rep.txt $TEMP_5 $PFX-mixcr-alig.vdjca
    mixcr assemble -f -t $NCORES -r $PFX-mixcr-asm.rep.txt $PFX-mixcr-alig.vdjca $PFX-mixcr-clone.clns

    rm $SEQ_1 $SEQ_2 $TEMP_3 $TEMP_3B $TEMP_3T $TEMP_3F $TEMP_4 $TEMP_5
    
fi


echo $(date) >> /job/docker_done
