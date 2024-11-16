#!/usr/bin/env bash
echo $(date) >> /job/docker_start

set -eu
shopt -s nullglob

#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate

source /job/config.txt

NCORES=$(nproc --all)
OUT=/output
PFX=/output/$ID

## Collect Input
echo $(date) "Collecting FASTQ files" >> /job/docker_script.log
IFS=';' read -ra FQ1SA <<< $FQ1S
FQ1SL=""
for i in "${FQ1SA[@]}"; do
    FID1=$(basename $i)
    FQ1SL+=" /input/$FID1"
done
IFS=';' read -ra FQ2SA <<< $FQ2S
FQ2SL=""
for i in "${FQ2SA[@]}"; do
    FID2=$(basename $i)
    FQ2SL+=" /input/$FID2"
done
tmp=$(tempfile -d /tmp) && rm $tmp
fq1cat=$tmp-cat_1.fq.gz
fq2cat=$tmp-cat_2.fq.gz
cat $FQ1SL > $fq1cat &
cat $FQ2SL > $fq2cat &
FQ1L=$fq1cat
FQ2L=$fq2cat

wait

## check and correct quality encoding
FORMAT=$(testformat.sh $FQ1L $FQ2L | cut -f1 | uniq)
if [ "$FORMAT" != "sanger" ]; then
    echo $(date) "Running Reformatting" >> /job/docker_script.log
    reformat.sh in=$FQ1L in2=$FQ2L out=/tmp/1.fq out2=/tmp/2.fq qin=auto qout=sanger t=$NCORES &>> /job/docker_script.log
    FQ1L=/tmp/1.fq
    FQ2L=/tmp/2.fq
fi

## trim adapters from reads
if [ "$ALIGN_TRIM" != "none" ]; then
    echo $(date) "Running Trimming" >> /job/docker_script.log
    bbduk.sh -Xmx1g t=$NCORES \
        in1=$FQ1L in2=$FQ2L out1=/tmp/1t.fq out2=/tmp/2t.fq ref=/opt/bbmap/resources/$ALIGN_TRIM.fa.gz \
        ktrim=r k=23 mink=11 hdist=1 ignorebadquality=t qin=33 \
        tpe tbo &>> /job/docker_script.log
    rm -f /tmp/1.fq /tmp/2.fq
    FQ1L=/tmp/1t.fq
    FQ2L=/tmp/2t.fq
fi

## umi consensus
echo $(date) "Running UMI Consensus" >> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon umi extract \
        $ALIGN_UMI \
        $FQ1L $FQ2L | \
$SENTIEON_INSTALL_DIR/bin/sentieon bwa mem \
        $ALIGN_BWA -p -C \
        -R "@RG\tID:$ID\tSM:$SAMPLE\tLB:$SAMPLE\tPL:ILLUMINA" \
        -t $NCORES $ALIGN_INDEX - 2>> /job/docker_script.log | \
$SENTIEON_INSTALL_DIR/bin/sentieon umi consensus \
        -o /tmp/12c.fq 2>> $PFX-umistat.txt
rm -f $fq1cat $fq2cat $FQ1L $FQ2L /tmp/1.fq /tmp/2.fq /tmp/1t.fq /tmp/2t.fq

## align and sort
echo $(date) "Running Alignment" >> /job/docker_script.log
( $SENTIEON_INSTALL_DIR/bin/bwa mem $ALIGN_BWA -p -C -R "@RG\tID:$ID\tSM:$SAMPLE\tLB:$SAMPLE\tPL:ILLUMINA" -t $NCORES $ALIGN_INDEX /tmp/12c.fq 2>> /job/docker_script.log || echo -n 'error' ) \
| $SENTIEON_INSTALL_DIR/bin/sentieon util sort -r $ALIGN_FASTA -o $PFX-raw.bam -t $NCORES --sam2bam --umi_post_process -i - &>> /job/docker_script.log
rm -f /tmp/12c.fq

## collect metrics
echo $(date) "Collecting QC Metrics" >> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $PFX-raw.bam \
				   --algo MeanQualityByCycle $PFX-meanqual.txt \
				   --algo QualDistribution   $PFX-qualdist.txt \
				   --algo GCBias --summary   $PFX-gcsummary.txt $PFX-gcbias.txt \
				   --algo AlignmentStat --adapter_seq '' $PFX-aligstat.txt \
				   --algo InsertSizeMetricAlgo $PFX-isize.txt &>> /job/docker_script.log

## indel realignment
echo $(date) "Indel Realignment" >> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $PFX-raw.bam --algo Realigner $ALIGN_INDEL $PFX-realigned.bam &>> /job/docker_script.log
rm $PFX-raw.bam*

## BQSR
echo $(date) "Computing BQSR tables" >> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $PFX-realigned.bam --algo QualCal $ALIGN_DBSNP $ALIGN_INDEL $PFX-recal_data.table &>> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $PFX-realigned.bam -q $PFX-recal_data.table --algo QualCal $ALIGN_DBSNP $ALIGN_INDEL $PFX-recal_data.table.post &>> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NCORES --algo QualCal --plot --before $PFX-recal_data.table --after $PFX-recal_data.table.post $PFX-recal.csv &>> /job/docker_script.log

## plots
echo $(date) "Plotting" >> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon plot bqsr -o $PFX-recal_plots.pdf $PFX-recal.csv &>> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon plot metrics -o $PFX-report.pdf gc=$PFX-gcbias.txt qd=$PFX-qualdist.txt mq=$PFX-meanqual.txt isize=$PFX-isize.txt &>> /job/docker_script.log

## genotyping 
echo $(date) "Generating GVCF" >> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $PFX-realigned.bam \
  --interval $ALIGN_GENOTYPE \
  --algo DNAscope --emit_mode gvcf \
  $PFX-genotype.gvcf &>> /job/docker_script.log

echo $(date) "Generating genotyping table" >> /job/docker_script.log
Rscript /code/pipe/common/genotype_calling.R \
  -i $PFX-genotype.gvcf \
  -o $PFX-genotype.csv &>> /job/docker_script.log

## clean-up
echo $(date) "Clean-up" >> /job/docker_script.log
mv $PFX-realigned.bam $PFX.bam
mv $PFX-realigned.bam.bai $PFX.bam.bai

echo $(date) "Finished." >> /job/docker_script.log
echo $(date) >> /job/docker_done
