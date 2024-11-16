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

PFX=/output/$ID

## collect input CRAMs/BAMs
echo $(date) "Collecting input files" >> /job/docker_script.log
BID=$(basename $PALN)
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

## Calculate QcMetrics (Sentieon)
if [ "$MISC_QCMETRICS" = true ]; then
    echo $(date) "Running QcMetrics" >> /job/docker_script.log
    if [ -z "$MISC_TARGETS" ]; then
        HSALGO=
    else
        HSALGO=-"-algo HsMetricAlgo --targets_list $MISC_TARGETS --baits_list $MISC_TARGETS $PFX-hsmetrics.txt"
    fi
    $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $BAM \
                                       $HSALGO \
                                       --algo SequenceArtifactMetricsAlgo $PFX-artifact \
                                       --algo MeanQualityByCycle $PFX-meanqual.txt \
                                       --algo QualDistribution $PFX-qualdist.txt \
                                       --algo GCBias --summary $PFX-gcsummary.txt $PFX-gcbias.txt \
                                       --algo AlignmentStat --adapter_seq '' $PFX-aligstat.txt \
                                       --algo InsertSizeMetricAlgo $PFX-isize.txt &>> /job/docker_script.log
fi

## Calculate CoverageMetrics (Sentieon)
if [ "$MISC_COVMETRICS" = true ]; then
    echo $(date) "Running CoverageMetrics" >> /job/docker_script.log
    if [ -z "$MISC_TARGETS" ]; then
        COVALGO="--algo WgsMetricsAlgo --include_unpaired true --min_map_qual 0 --min_base_qual 0 --coverage_cap 300 $PFX-wgsmetrics.txt"
    else
        if [ "$MISC_COVPERBASE" =  true ]; then
            COVARGS=--omit_base_output
        else
            COVARGS=
        fi
        COVALGO="--interval $MISC_TARGETS --algo CoverageMetrics $COVARGS $PFX-covmetrics"
    fi
    $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $BAM \
                                       $COVALGO &>> /job/docker_script.log
    if [[ -s $PFX-covmetrics ]]; then
        echo $(date) "Compressing covmetrics" >> /job/docker_script.log
        pigz -p $NCORES $PFX-covmetrics
    fi
fi

if [ "$MISC_GENOTYPE" = true ]; then
    ## genotyping
    echo $(date) "Generating GVCF" >> /job/docker_script.log
    $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $BAM \
                                       --interval $ALIGN_GENOTYPE \
                                       --algo DNAscope --emit_mode gvcf \
                                       $PFX-genotype.gvcf &>> /job/docker_script.log

    echo $(date) "Generating genotyping table" >> /job/docker_script.log
    /usr/bin/Rscript /code/pipe/common/genotype_calling.R \
            -i $PFX-genotype.gvcf \
            -o $PFX-genotype.csv &>> /job/docker_script.log
fi

## Run Virmer
if [ "$MISC_VIRMER" = true ]; then
    # extract unmapped reads
    samtools view -@4 -u -f 12 -F 256 $BAM > $PFX-tmp.bam
    # convert bam to fastq
    bamtofastq \
        collate=1 \
        exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
        filename=$PFX-tmp.bam \
        gz=0 \
        inputformat=bam \
        F=$PFX-tmp_1.fq \
        F2=$PFX-tmp_2.fq \
        S=$PFX-tmp_s.fq \
        O=$PFX-tmp_o1.fq \
        O2=$PFX-tmp_o2.fq \
        tryoq=1 &>> /job/docker_script.log
    # quantify virus kmers
    bbduk.sh -Xmx8096M in=$PFX-tmp_1.fq in2=$PFX-tmp_2.fq ref=$MISC_VIRUSFASTA \
             $MISC_VIRMERARGS \
             stats=$PFX-virus.txt &>> /job/docker_script.log
    # cleanup
    rm $PFX-tmp*
fi

if [ "$MISC_MERGE" = true ]; then
    echo $(date) "Merging QC results" >> /job/docker_script.log
    /usr/bin/Rscript /code/pipe/cords-misc/merge_qc.R $BID $ID
fi

echo $(date) >> /job/docker_done
