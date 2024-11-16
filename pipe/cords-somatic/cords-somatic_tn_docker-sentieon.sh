## Run Co-realignment
if [ "$SOMATIC_REALIGN" = true ]; then
    echo $(date) "Running Co-realignment" >> /job/docker_script.log
    $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES $IS \
                                       --algo Realigner $ALIGN_INDEL $PFX-corealigned.bam &>> /job/docker_script.log
    IS="-i $PFX-corealigned.bam"
fi

## Run TNscope (improved GATK3 MuTect2)
if [ "$SOMATIC_TNSCOPE" = true ]; then
    echo $(date) "Running TNscope" >> /job/docker_script.log
    $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA \
        -t $NCORES $IS \
        --algo TNscope \
        --disable_detector sv \
        $SOMATIC_TNSCOPEARGS \
        --dbsnp $SOMATIC_DBSNP \
        --tumor_sample $BAMT_SAMPLE --normal_sample $BAMN_SAMPLE \
        $PFX-somatic-tnscope.vcf.gz &>> /job/docker_script.log
fi

## Run TNscope (in tumor-in-normal mode)
if [ "$SOMATIC_TINSCOPE" = true ]; then
    echo $(date) "Running TINscope" >> /job/docker_script.log
    $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA \
        -t $NCORES $IS \
        --algo TNscope \
        --disable_detector sv \
        $SOMATIC_TINSCOPEARGS \
        --dbsnp $SOMATIC_DBSNP \
        --tumor_sample $BAMT_SAMPLE --normal_sample $BAMN_SAMPLE \
        $PFX-somatic-tinscope.vcf.gz &>> /job/docker_script.log
fi

## Run DNAscope (improved GATK4 HaplotypeCaller)
if [ "$SOMATIC_DNASCOPE" = true ]; then
    echo $(date) "Running DNAscope" >> /job/docker_script.log
    $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA \
        -t $NCORES $IS \
        --algo DNAscope \
        --annotation QD --annotation MQ --annotation MQRankSum --annotation FS --annotation SOR --annotation MQ0 \
        $SOMATIC_DNASCOPEARGS \
        -d $SOMATIC_DBSNP \
        $PFX-somatic-dnascope.vcf.gz &>> /job/docker_script.log
fi

# clean-up
rm -f $PFX-corealigned.bam*
rm -f /output/tmp*
rm -rf /output/strelka-somatic
