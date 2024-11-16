docker_script_log=/job/docker_script.log
MEMava=$(($MEM*1000*95/1000)) # m
WHEREISJAR=`gatk --list 2>&1 >/dev/null | head -1|tr ' ' '\n'|tail -1`
# save 5% of memory for operating system, and change to M
MEMsp=$(($MEM * 95 * 1000 / 100 / $NCORES))

# collect input CRAMs/BAMs
echo $(date) "Collecting input files" >> $docker_script_log
IS=""
IFS=';' read -ra ALNA <<< $ALN
for i in "${ALNA[@]}"; do
    BID=$(basename $i)
    if [ -f "/input/$BID/$BID.cram" ]; then
        BAM="/input/$BID/$BID.cram"
        if [[ ! -f "/input/$BID/$BID.crai" && ! -f "/input/$BID/$BID.cram.crai" ]]; then
          samtools index -@ $NCORES $BAM /input/$BID/$BID.cram.crai
        fi
    elif [ -f "/input/$BID/$BID.bam" ]; then
        BAM="/input/$BID/$BID.bam"
        if [[ ! -f "/input/$BID/$BID.bai" && ! -f "/input/$BID/$BID.bam.bai" ]]; then
          samtools index -@ $NCORES $BAM /input/$BID/$BID.bam.bai
        fi
    fi
    # fixrg is not available for open-source version
    if [ "$ALIGN_FIXRG" = true ]; then
        /code/pipe/common/fixrg.sh "$BID-FIX" $ID $BAM
    fi

    IS+=" -I $BAM"
done

# mark dupcalites
echo $(date) "Removing Duplicates" >> $docker_script_log
gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" MarkDuplicates \
     $IS \
     -O $PFX-dedup.bam \
     -M $PFX-dedup_metrics.txt \
     --REMOVE_DUPLICATES \
     --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
     -R $ALIGN_FASTA \
     >> $docker_script_log 2>&1

# BQSR
echo $(date) "BQSR - calculation" >> $docker_script_log
gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" BaseRecalibratorSpark \
    -R $ALIGN_FASTA \
    -I $PFX-dedup.bam \
    $ALIGN_DBSNP \
    $ALIGN_INDEL \
    -O $PFX-recal_data.table \
    --use-original-qualities \
    >> $docker_script_log 2>&1

# applyBQSR
echo $(date) "BQSR - apply recalibration" >> $docker_script_log
gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" ApplyBQSRSpark \
    -R $ALIGN_FASTA \
    -I $PFX-dedup.bam \
    -bqsr $PFX-recal_data.table \
    -O $PFX-recaled.bam \
    --use-original-qualities \
    --create-output-bam-index \
    >> $docker_script_log 2>&1

## clean-up
echo $(date) "Clean-up" >> $docker_script_log
mv $PFX-recaled.bam $PFX.bam
if [[ -f "$PFX-recaled.bam.bai" ]]; then
    mv $PFX-recaled.bam.bai $PFX.bam.bai
fi
if [[ -f "$PFX-recaled.bai" ]]; then
    mv $PFX-recaled.bai $PFX.bai
fi
if [[ -f "$PFX-recaled.bam.sbi" ]]; then
    mv $PFX-recaled.bam.sbi $PFX.bam.sbi
fi
rm $PFX-dedup.bam*
