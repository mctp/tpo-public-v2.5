docker_script_log=/job/docker_script.log
MEMava=$(($MEM*1000*95/1000)) # m
WHEREISJAR=`gatk --list 2>&1 >/dev/null | head -1|tr ' ' '\n'|tail -1`

# save 5% of memory for operating system, and change to M
MEMsp=$(($MEM * 95 * 1000 / 100 / $NCORES))

## Run GATK4 MuTect2 (tumor-in-normal mode)
if [ "$SOMATIC_TINSCOPE" = true ] ; then
    echo "WARNING: MuTect2 (tumor-in-normal mode) is not available for open-source version. Even though TINSCOPE is set to be true in the config file, nothing will be output in this mode." >> $docker_script_log
fi

# create inteval_list
# ScatterIntervalsByNs
echo $(date) "ScatterIntervalsByNs" >> $docker_script_log
gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" ScatterIntervalsByNs \
    -R $ALIGN_FASTA \
    -O /output/byN.interval_list \
    --OUTPUT_TYPE ACGT \
    &>> $docker_script_log

# remove reference besides chromosome
grep -w -e 'chr1\|chr2\|chr3\|chr4\|chr5\|chr6\|chr7\|chr8\|chr9\|chr10\|chr11\|chr12\|chr13\|chr14\|chr15\|chr16\|chr17\|chr18\|chr19\|chr20\|chr21\|chr22\|chrX\|chrY' /output/byN.interval_list > /output/byN.interval_list.tmp
mv /output/byN.interval_list.tmp /output/byN.interval_list

# SplitIntervals
INTERVAL=/output/intervals
mkdir $INTERVAL
echo $(date) "SplitIntervals" >> $docker_script_log
gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" SplitIntervals \
    -R $ALIGN_FASTA \
    -L /output/byN.interval_list \
    -scatter-count $NCORES \
    -O $INTERVAL \
    &>> $docker_script_log

## Run GATK4 MuTect2
if [ "$SOMATIC_TNSCOPE" = true ] ; then
    echo $(date) "Running MuTect2" >> $docker_script_log
    MUT_TMP=/output/MUT_TMP
    mkdir $MUT_TMP
    ## parallel running of Mutect2
    parallel "-j${NCORES}" \
    java -Xmx${MEMsp}m -XX:ConcGCThreads=1 -jar $WHEREISJAR Mutect2 \
        -R $ALIGN_FASTA \
        $IS \
        $SOMATIC_TNSCOPEARGS \
        --germline-resource $SOMATIC_GERMLINEDB \
        -normal $BAMN_SAMPLE \
        -O $MUT_TMP/$ID-{/.}-unfiltered.vcf.gz \
        --f1r2-tar-gz $MUT_TMP/$ID-{/.}.f1r2.tar.gz \
        -L {} \
        --native-pair-hmm-threads 1 \
        "&>>" $MUT_TMP/parallel-mutect-{/.}.log \
        ::: $INTERVAL/*.interval_list

    cat $MUT_TMP/parallel-mutect-*.log >> $docker_script_log
    rm $MUT_TMP/parallel-mutect-*.log

    ## get pileups and calculate the contaminations
    ### for tumor
    echo $(date) "Get pileups for tumor" >> $docker_script_log
    parallel "-j${NCORES}" \
    java -Xmx${MEMsp}m -XX:ConcGCThreads=1 -jar $WHEREISJAR GetPileupSummaries \
        -R $ALIGN_FASTA \
        -I $BAMT \
        -V $SOMATIC_GERMLINEDB \
        -L {} \
        -O $MUT_TMP/$ID-{/.}-tumor-pileups.table \
        "&>>" $MUT_TMP/parallel-pileup-{/.}.log \
        ::: $INTERVAL/*.interval_list

    cat $MUT_TMP/parallel-pileup-*.log >> $docker_script_log
    rm $MUT_TMP/parallel-pileup-*.log

    ## get pileups and calculate the contaminations
    ### for normal
    echo $(date) "Get pileups for normal" >> $docker_script_log
    parallel "-j${NCORES}" \
    java -Xmx${MEMsp}m -XX:ConcGCThreads=1 -jar $WHEREISJAR GetPileupSummaries \
        -R $ALIGN_FASTA \
        -I $BAMN \
        -V $SOMATIC_GERMLINEDB \
        -L {} \
        -O $MUT_TMP/$ID-{/.}-normal-pileups.table \
        "&>>" $MUT_TMP/parallel-pileup-{/.}.log \
        ::: $INTERVAL/*.interval_list

    cat $MUT_TMP/parallel-pileup-*.log >> $docker_script_log
    rm $MUT_TMP/parallel-pileup-*.log

    ## learn read orientation model
    echo $(date) "Learn read orientation models" >> $docker_script_log
    F1R2S=$(ls -1 $MUT_TMP/$ID-*.f1r2.tar.gz | sed "s#^#-I #" |tr '\n' ' ')
    gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" LearnReadOrientationModel \
        $F1R2S \
        -O $MUT_TMP/$ID-priors.tar.gz \
        &>> $docker_script_log

    ## Merge unfiltered VCFs
    echo $(date) "Merge unfiltered VCFs" >> $docker_script_log
    UNFILTERED_VCFS=$(ls -1 $MUT_TMP/$ID-*-unfiltered.vcf.gz | sed "s#^#-I #" |tr '\n' ' ')
    gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" MergeVcfs \
        $UNFILTERED_VCFS \
        -O $MUT_TMP/$ID-unfiltered.vcf.gz \
        &>> $docker_script_log

    ## MergeStats
    echo $(date) "Merge Mutect2 stats" >> $docker_script_log
    MUTECT2STATS=$(ls -1 $MUT_TMP/$ID-*-unfiltered.vcf.gz.stats | sed "s#^#-stats #" |tr '\n' ' ')
    gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" MergeMutectStats \
        $MUTECT2STATS \
        -O $MUT_TMP/$ID-unfiltered.vcf.gz.stats \
        &>> $docker_script_log

    ## MergePileupSummaries
    ## for tumor
    echo $(date) "Merge pipeup summaries for tumor" >> $docker_script_log
    PILEUPS=$(ls -1 $MUT_TMP/$ID-*-tumor-pileups.table | sed "s#^#-I #" |tr '\n' ' ')
    gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" GatherPileupSummaries \
        $PILEUPS \
        --sequence-dictionary $(echo $ALIGN_FASTA | sed 's/\.[^.]*$/.dict/') \
        -O $MUT_TMP/$ID-tumor-pileups.table \
        &>> $docker_script_log

    ## MergePileupSummaries
    ## for normal
    echo $(date) "Merge pipeup summaries for normal" >> $docker_script_log
    PILEUPS=$(ls -1 $MUT_TMP/$ID-*-normal-pileups.table | sed "s#^#-I #" |tr '\n' ' ')
    gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" GatherPileupSummaries \
        $PILEUPS \
        --sequence-dictionary $(echo $ALIGN_FASTA | sed 's/\.[^.]*$/.dict/') \
        -O $MUT_TMP/$ID-normal-pileups.table \
        &>> $docker_script_log

    ## Calculate contaimination
    echo $(date) "Calculate contamination" >> $docker_script_log
    gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" CalculateContamination \
        -I $MUT_TMP/$ID-tumor-pileups.table \
        -matched $MUT_TMP/$ID-normal-pileups.table \
        --tumor-segmentation $MUT_TMP/$ID-segments.table \
        -O $MUT_TMP/$ID-contamination.table \
        &>> $docker_script_log

    ## filtering
    echo $(date) "Filter mutect calls" >> $docker_script_log
    gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" FilterMutectCalls \
        -V $MUT_TMP/$ID-unfiltered.vcf.gz \
        -R $ALIGN_FASTA \
        -O $PFX-filtered.vcf.gz \
        --contamination-table $MUT_TMP/$ID-contamination.table \
        --tumor-segmentation $MUT_TMP/$ID-segments.table \
        -ob-priors $MUT_TMP/$ID-priors.tar.gz \
        --stats $MUT_TMP/$ID-unfiltered.vcf.gz.stats \
        --filtering-stats $PFX-filtered.vcf.gz.stats \
        &>> $docker_script_log

    mv $PFX-filtered.vcf.gz.stats $PFX-somatic-mutect2.vcf.gz.stats
    mv $PFX-filtered.vcf.gz $PFX-somatic-mutect2.vcf.gz
    mv $PFX-filtered.vcf.gz.tbi $PFX-somatic-mutect2.vcf.gz.tbi

    ## remove tmp files
    rm -r $MUT_TMP
fi

## Run GATK4 HaplotypeCaller
if [ "$SOMATIC_DNASCOPE" = true ] ; then
    echo $(date) "Running haplotypeCaller" >> $docker_script_log
    HAP_TMP=/output/HAP_TMP
    mkdir $HAP_TMP
    ## parallel running of HaplotypeCaller
    parallel "-j${NCORES}" \
    java -Xmx${MEMsp}m -XX:ConcGCThreads=1 -jar $WHEREISJAR HaplotypeCaller \
        -R $ALIGN_FASTA \
        $IS \
        --annotation QualByDepth --annotation MappingQuality --annotation MappingQualityRankSumTest --annotation FisherStrand --annotation StrandOddsRatio --annotation MappingQualityZero \
        $SOMATIC_DNASCOPEARGS \
        -D $SOMATIC_DBSNP \
        -O $HAP_TMP/$ID-{/.}-somatic-haplotypecaller.vcf.gz \
        -L {} \
        --native-pair-hmm-threads 1 \
        "&>>" $HAP_TMP/parallel-{/.}.log \
        ::: $INTERVAL/*.interval_list

    cat $HAP_TMP/parallel-*.log >> $docker_script_log

    ## store the path of all vcf.gz
    ls -1 $HAP_TMP/*.vcf.gz > $HAP_TMP/vcf.list
    gatk --java-options "-Xmx${MEMava}m -XX:ConcGCThreads=${NCORES}" MergeVcfs \
        -I $HAP_TMP/vcf.list \
        -O $PFX-somatic-haplotypecaller.vcf.gz

    ## remove tmp files
    rm -r $HAP_TMP
fi

## remove interval tmp files
rm /output/byN.interval_list
rm -r $INTERVAL
