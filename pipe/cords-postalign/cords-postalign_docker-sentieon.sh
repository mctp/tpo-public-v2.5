## collect input CRAMs/BAMs
echo $(date) "Collecting input files" >> /job/docker_script.log
IS=""
QS=""
IFS=';' read -ra ALNA <<< $ALN
for i in "${ALNA[@]}"; do
    BID=$(basename $i)
    if [ -f "/input/$BID/$BID.cram" ]; then
        BAM="/input/$BID/$BID.cram"
        if [ ! -f "/input/$BID/$BID.crai" -a ! -f "/input/$BID/$BID.cram.crai" ]; then
          samtools index -@ $NCORES $BAM /input/$BID/$BID.cram.crai
        fi
    elif [ -f "/input/$BID/$BID.bam" ]; then
        BAM="/input/$BID/$BID.bam"
        if [ ! -f "/input/$BID/$BID.bai" -a ! -f "/input/$BID/$BID.bam.bai" ]; then
          samtools index -@ $NCORES $BAM /input/$BID/$BID.bam.bai
        fi
    fi
    if [ "$ALIGN_FIXRG" = true ]; then
        /code/pipe/common/fixrg.sh "$BID-FIX" $ID $BAM
    fi
    IS+=" -i $BAM"
    if [ -f "/input/$BID/$BID-recal_data.table" ]; then
        QS+=" -q /input/$BID/$BID-recal_data.table"
    fi
done

echo $(date) "Removing Duplicates" >> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES $IS $QS \
                --algo LocusCollector --fun score_info $PFX-dup-score.txt &>> /job/docker_script.log

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES $IS $QS \
                --algo Dedup --optical_dup_pix_dist 2500 --rmdup --score_info $PFX-dup-score.txt --metrics $PFX-dedup_metrics.txt $PFX.bam &>> /job/docker_script.log

## clean-up
rm -f $PFX-dup-score.txt
rm -f /output/*.idx
