docker_script_log=/job/docker_script.log
MEMava=$(($MEM*1000*95/1000)) # m
WHEREISJAR=`gatk --list 2>&1 >/dev/null | head -1|tr ' ' '\n'|tail -1`
# save 5% of memory for operating system, and change to M
MEMsp=$(($MEM * 95 * 1000 / 100 / $NCORES))

# align and sort
echo $(date) "Running Alignment" >> $docker_script_log
(bwa mem \
    $ALIGN_BWA \
    -R "@RG\tID:$ID\tSM:$SAMPLE\tLB:$SAMPLE\tPL:ILLUMINA" \
    -t $NCORES \
    $ALIGN_INDEX \
    $FQ1L \
    $FQ2L \
    -o /dev/stdout \
    2>> $docker_script_log || echo -n 'error') \
| samtools sort \
    --reference $ALIGN_FASTA \
    -o $PFX-sorted.cram##idx##$PFX-sorted.cram.crai \
    -@ $NCORES \
    -l 6 \
    --write-index \
    - \
    >> $docker_script_log 2>&1
rm -f /tmp/1.fq /tmp/2.fq /tmp/1t.fq /tmp/2t.fq /tmp/tmp_*.config

# collecting QC metrics
echo $(date) "Collecting QC Metrics" >> $docker_script_log
# Spark running and &+wait parallel running
echo $(date) "MeanQualityByCycleSpark" >> /job/metrics1.log
gatk --java-options "-Xmx$(($MEMava/5))m -XX:ConcGCThreads=$(($NCORES/5))" MeanQualityByCycleSpark \
    -R $ALIGN_FASTA \
    -I $PFX-sorted.cram \
    -O $PFX-meanqual.txt \
    --chart $PFX-meanqual.pdf \
    >> /job/metrics1.log 2>&1 &

echo $(date) "QualityScoreDistributionSpark" >> /job/metrics2.log
gatk --java-options "-Xmx$(($MEMava/5))m -XX:ConcGCThreads=$(($NCORES/5))" QualityScoreDistributionSpark \
    -R $ALIGN_FASTA \
    -I $PFX-sorted.cram \
    -O $PFX-qualdist.txt \
    --chart $PFX-qualdist.pdf \
    >> /job/metrics2.log 2>&1 &

echo $(date) "CollectGcBiasMetrics" >> /job/metrics3.log
gatk --java-options "-Xmx$(($MEMava/5))m -XX:ConcGCThreads=$(($NCORES/5))" CollectGcBiasMetrics \
    -R $ALIGN_FASTA \
    -I $PFX-sorted.cram \
    -O $PFX-gcbias.txt \
    -CHART $PFX-gcbias.pdf \
    -S $PFX-gcsummary.txt \
    >> /job/metrics3.log 2>&1 &

echo $(date) "CollectAlignmentSummaryMetrics" >> /job/metrics4.log
gatk --java-options "-Xmx$(($MEMava/5))m -XX:ConcGCThreads=$(($NCORES/5))" CollectAlignmentSummaryMetrics \
    -R $ALIGN_FASTA \
    -I $PFX-sorted.cram \
    -O $PFX-aligstat.txt \
    --ADAPTER_SEQUENCE '' \
    >> /job/metrics4.log 2>&1 &

echo $(date) "CollectInsertSizeMetrics" >> /job/metrics5.log
gatk --java-options "-Xmx$(($MEMava/5))m -XX:ConcGCThreads=$((NCORES/5))" CollectInsertSizeMetricsSpark \
    -R $ALIGN_FASTA \
    -I $PFX-sorted.cram \
    -O $PFX-isize.txt \
    -H $PFX-isize.pdf \
    >> /job/metrics5.log 2>&1 &
wait
cat /job/metrics1.log /job/metrics2.log /job/metrics3.log /job/metrics4.log /job/metrics5.log >> $docker_script_log
rm /job/metrics1.log /job/metrics2.log /job/metrics3.log /job/metrics4.log /job/metrics5.log

# genotyping
echo $(date) "Generating GVCF" >> $docker_script_log
gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" HaplotypeCaller \
    -R $ALIGN_FASTA \
    -I $PFX-sorted.cram \
    -L $ALIGN_GENOTYPE \
    -ERC GVCF \
    -O $PFX-genotype.gvcf \
    >> $docker_script_log 2>&1

echo $(date) "Generating genotyping table" >> $docker_script_log
/usr/bin/Rscript /code/pipe/common/genotype_calling.R \
    -i $PFX-genotype.gvcf \
    -o $PFX-genotype.csv \
    >> $docker_script_log 2>&1

# clean
echo $(date) "Clean-up" >> $docker_script_log
mv $PFX-sorted.cram $PFX.cram
mv $PFX-sorted.cram.crai $PFX.cram.crai
