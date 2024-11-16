#!/usr/bin/env bash
echo $(date) >> /job/docker_start

set -eu
shopt -s nullglob


DNASCOPE_CALL () {
    $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA $1 \
        -t $NCORES $IS \
        --algo DNAscope \
        --annotation QD --annotation MQ --annotation MQRankSum --annotation FS --annotation SOR --annotation MQ0 \
        $GERMLINE_DNASCOPEARGS \
        $2 \
        -d $GERMLINE_DBSNP \
        $PFX-$3-germline-dnascope.vcf.gz &>> /job/docker_script.log
}

#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate

source /job/config.txt
NCORES=$(nproc --all)

PFX=/output/$ID

## collect input CRAMs/BAMs
echo $(date) "Collecting input files" >> /job/docker_script.log
BID=$(basename $ALN)
if [ -f "/input/$BID/$BID.cram" ]; then
    BAM="/input/$BID/$BID.cram"
elif [ -f "/input/$BID/$BID.bam" ]; then
    BAM="/input/$BID/$BID.bam"
fi
if [ ! -f "/input/$BID/$BID.cram.crai" ] && [ ! -f "/input/$BID/$BID.bam.bai" ] && \
   [ ! -f "/input/$BID/$BID.bai" ] && [ ! -f "/input/$BID/$BID.crai" ]; then
    $SENTIEON_INSTALL_DIR/bin/sentieon util index $BAM
fi
IS="-i $BAM"

## Determining Sex
if [ -z "$SEX" ]; then
    if [ -z "$GERMLINE_SEXSNP" ]; then
        echo $(date) "Sexing based on SNPs" >> /job/docker_script.log
        $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $BAM \
        --interval $GERMLINE_SEXSNP \
        --algo DNAscope --emit_mode gvcf \
        $PFX-sex.gvcf &>> /job/docker_script.log
        /usr/bin/Rscript /code/pipe/common/sex_calling.R \
            -s $PFX-sex.gvcf \
            -o $PFX-sex.txt &>> /job/docker_script.log
    else
        /usr/bin/Rscript /code/pipe/common/sex_calling.R \
            -b $BAM \
            -o $PFX-sex.txt &>> /job/docker_script.log
    fi
    source $PFX-sex.txt
else
    echo "SEX=$SEX" > $PFX-sex.txt
fi

## Run Co-realignment
if [ "$GERMLINE_REALIGN" = true ]; then
    echo $(date) "Running Realignment" >> /job/docker_script.log
    $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES $IS \
	                               --algo Realigner $ALIGN_INDEL $PFX-realigned.bam &>> /job/docker_script.log
    IS="-i $PFX-realigned.bam"
fi

## Run DNAscope (improved GATK4 HaplotypeCaller)
if [ "$GERMLINE_DNASCOPE" = true ] ; then
    echo "Running DNAscope" >> /job/docker_script.log
    if [ "$SEX" = "XX" ]; then
        DNASCOPE_CALL "--interval $GERMLINE_INTERVALS-autosome.bed" "--ploidy 2" autosome
        DNASCOPE_CALL "--interval $GERMLINE_INTERVALS-x-par1.bed" "--ploidy 2" x-par1
        DNASCOPE_CALL "--interval $GERMLINE_INTERVALS-x-nonpar.bed" "--ploidy 2" x-nonpar
        DNASCOPE_CALL "--interval $GERMLINE_INTERVALS-x-par2.bed" "--ploidy 2" x-par2
    fi
    if [ "$SEX" = "XY" ]; then
        DNASCOPE_CALL "--interval $GERMLINE_INTERVALS-autosome.bed" "--ploidy 2" autosome
        DNASCOPE_CALL "--interval $GERMLINE_INTERVALS-x-par1.bed" "--ploidy 2" x-par1
        DNASCOPE_CALL "--interval $GERMLINE_INTERVALS-x-nonpar.bed" "--ploidy 1" x-nonpar
        DNASCOPE_CALL "--interval $GERMLINE_INTERVALS-x-par2.bed" "--ploidy 2" x-par2
    fi
    ## force is used because dnascope does not provide FORMAT/PGT FORMAT/PID for ploidy=1
    bcftools concat --naive-force \
             $PFX-autosome-germline-dnascope.vcf.gz \
             $PFX-x-par1-germline-dnascope.vcf.gz \
             $PFX-x-nonpar-germline-dnascope.vcf.gz \
             $PFX-x-par2-germline-dnascope.vcf.gz \
             -o $PFX-germline-dnascope.vcf.gz &>> /job/docker_script.log
    tabix $PFX-germline-dnascope.vcf.gz
fi

## clean-up
rm -f $PFX-realigned.bam*
rm -f $PFX-{autosome-,x-}*
rm -f $PFX-sex.gvcf $PFX-sex.gvcf.idx

echo $(date) >> /job/docker_done
