## align and sort
echo $(date) "Running Alignment" >> /job/docker_script.log
( $SENTIEON_INSTALL_DIR/bin/bwa mem $ALIGN_BWA -R "@RG\tID:$ID\tSM:$SAMPLE\tLB:$SAMPLE\tPL:ILLUMINA" -t $NCORES $ALIGN_INDEX $FQ1L $FQ2L 2>> /job/docker_script.log || echo -n 'error' ) \
| $SENTIEON_INSTALL_DIR/bin/sentieon util sort -r $ALIGN_FASTA -o $PFX-raw.bam -t $NCORES --sam2bam -i - &>> /job/docker_script.log
rm -f /tmp/1.fq /tmp/2.fq /tmp/1t.fq /tmp/2t.fq

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
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $PFX-raw.bam --algo Realigner $ALIGN_INDEL --cram_write_options version=3.0,compressor=rans $PFX-realigned.cram &>> /job/docker_script.log

## BQSR
echo $(date) "Computing BQSR tables" >> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $PFX-realigned.cram --algo QualCal $ALIGN_DBSNP $ALIGN_INDEL $PFX-recal_data.table &>> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $PFX-realigned.cram -q $PFX-recal_data.table --algo QualCal $ALIGN_DBSNP $ALIGN_INDEL $PFX-recal_data.table.post &>> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $NCORES --algo QualCal --plot --before $PFX-recal_data.table --after $PFX-recal_data.table.post $PFX-recal.csv &>> /job/docker_script.log

## plots
echo $(date) "Plotting" >> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon plot bqsr -o $PFX-recal_plots.pdf $PFX-recal.csv &>> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon plot metrics -o $PFX-report.pdf gc=$PFX-gcbias.txt qd=$PFX-qualdist.txt mq=$PFX-meanqual.txt isize=$PFX-isize.txt &>> /job/docker_script.log

## genotyping 
echo $(date) "Generating GVCF" >> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $PFX-realigned.cram \
  --interval $ALIGN_GENOTYPE \
  --algo DNAscope --emit_mode gvcf \
  $PFX-genotype.gvcf &>> /job/docker_script.log

echo $(date) "Generating genotyping table" >> /job/docker_script.log
/usr/bin/Rscript /code/pipe/common/genotype_calling.R \
    -i $PFX-genotype.gvcf \
    -o $PFX-genotype.csv &>> /job/docker_script.log
 
## clean-up
echo $(date) "Clean-up" >> /job/docker_script.log
mv $PFX-realigned.cram $PFX.cram
mv $PFX-realigned.cram.crai $PFX.cram.crai
mv $PFX-realigned.cram.bai $PFX.cram.bai
rm $PFX-raw.bam*

