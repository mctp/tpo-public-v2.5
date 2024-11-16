#!/usr/bin/env bash
echo $(date) >> /job/docker_start

set -eu
shopt -s nullglob

#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate crisp

source /job/config.txt
NCORES=$(nproc --all)

PFX=/output/$ID

## collect input CRAMs/BAMs
echo $(date) "Collecting input files" >> /job/docker_script.log
BID=$(basename $ALNT)
if [ -f "/input/$BID/$BID.cram" ]; then
    BAMT="/input/$BID/$BID.cram"
elif [ -f "/input/$BID/$BID.bam" ]; then
    BAMT="/input/$BID/$BID.bam"
fi
if [ ! -f "/input/$BID/$BID.cram.crai" ] && [ ! -f "/input/$BID/$BID.bam.bai" ] && \
   [ ! -f "/input/$BID/$BID.bai" ] && [ ! -f "/input/$BID/$BID.crai" ]; then
    $SENTIEON_INSTALL_DIR/bin/sentieon util index $BAMT
fi
IS="-i $BAMT"

## Remove Duplicates
echo $(date) "Removing Duplicates" >> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES $IS \
                                   --algo LocusCollector --fun score_info $PFX-dup-score.txt &>> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES $IS \
                                   --algo Dedup --optical_dup_pix_dist 2500 --rmdup --score_info $PFX-dup-score.txt --metrics $PFX-dedup_metrics.txt $PFX-dedup.bam &>> /job/docker_script.log

## Split Reads
echo $(date) "Splitting Reads" >> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $PFX-dedup.bam --algo RNASplitReadsAtJunction --reassign_mapq 255:60 $PFX-split.bam &>> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $PFX-split.bam --algo QualCal $ALIGN_DBSNP $ALIGN_INDEL $PFX-recal_data.table &>> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $PFX-split.bam -q $PFX-recal_data.table --algo QualCal $ALIGN_DBSNP $ALIGN_INDEL $PFX-recal_data.table.post &>> /job/docker_script.log

## Run DNAscope (improved GATK4 HaplotypeCaller)
if [ "$GERMLINE_DNASCOPE" = true ] ; then
    echo $(date) "Calling Variants" >> /job/docker_script.log
    $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA \
        -t $NCORES -i $PFX-split.bam -q $PFX-recal_data.table \
        --algo DNAscope \
        $GERMLINE_DNASCOPEARGS \
        -d $GERMLINE_DBSNP \
        $PFX-germline-dnascope.vcf.gz &>> /job/docker_script.log
fi

## clean-up
rm -f $PFX-recal_data.table $PFX-recal_data.table.post
rm -f $PFX-dedup.bam* $PFX-split.bam*
rm -f $PFX-dup-score.txt
rm -f /output/*.idx

echo $(date) >> /job/docker_done
