#!/usr/bin/env bash
echo $(date) >> /job/docker_start

set -eu
shopt -s nullglob

#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate structural

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
BID=$(basename $ALNN)
if [ -f "/input/$BID/$BID.cram" ]; then
    BAMN="/input/$BID/$BID.cram"
elif [ -f "/input/$BID/$BID.bam" ]; then
    BAMN="/input/$BID/$BID.bam"
fi
if [ ! -f "/input/$BID/$BID.cram.crai" ] && [ ! -f "/input/$BID/$BID.bam.bai" ] && \
   [ ! -f "/input/$BID/$BID.bai" ] && [ ! -f "/input/$BID/$BID.crai" ]; then
    $SENTIEON_INSTALL_DIR/bin/sentieon util index $BAMN
fi
IS="-i $BAMT -i $BAMN"
BAMT_SAMPLE=$(samtools view -H $BAMT | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
BAMN_SAMPLE=$(samtools view -H $BAMN | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)

## Run TNScope
if [ "$STRUCTURAL_TNSCOPE" = true ]; then
    echo $(date) "Running TNscope" >> /job/docker_script.log
    $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA \
        -t $NCORES $IS \
        --algo TNscope \
        --disable_detector snv_indel \
        $STRUCTURAL_TNSCOPEARGS \
        --tumor_sample $BAMT_SAMPLE --normal_sample $BAMN_SAMPLE \
        $PFX-structural-tnscope.vcf.gz &>> /job/docker_script.log
fi

## Run Manta
if [ "$STRUCTURAL_MANTA" = true ]; then
    echo $(date) "Running Manta Somatic" >> /job/docker_script.log
    configManta.py \
        --tumorBam $BAMT \
        --normalBam $BAMN \
        --referenceFasta $ALIGN_FASTA \
        --runDir /output/manta-somatic \
        $STRUCTURAL_MANTAARGS
    /output/manta-somatic/runWorkflow.py -j $NCORES &>> /job/docker_script.log

    for fn in /output/manta-somatic/results/variants/*
    do
        lcfn=$(echo $(basename $fn) | tr '[:upper:]' '[:lower:]')
        mv "$fn" $PFX-structural-manta-$lcfn
    done
    
    ## clean-up
    rm -rf /output/manta-somatic
fi

## Run Manta
if [ "$STRUCTURAL_SVABA" = true ]; then
    echo $(date) "Running SvABA Somatic" >> /job/docker_script.log

    mkdir -p /output/svaba-tmp
    ln -sfn $ALIGN_FASTA /output/svaba-tmp/reference.fa
    ln -sfn $ALIGN_FASTA.fai /output/svaba-tmp/reference.fa.fai
    ln -sfn $ALIGN_INDEX.dict /output/svaba-tmp/reference.fa.dict
    ln -sfn $ALIGN_INDEX.alt /output/svaba-tmp/reference.fa.alt
    ln -sfn $ALIGN_INDEX.ann /output/svaba-tmp/reference.fa.ann
    ln -sfn $ALIGN_INDEX.pac /output/svaba-tmp/reference.fa.pac
    ln -sfn $ALIGN_INDEX.amb /output/svaba-tmp/reference.fa.amb
    ln -sfn $ALIGN_INDEX.bwt /output/svaba-tmp/reference.fa.bwt
    ln -sfn $ALIGN_INDEX.sa /output/svaba-tmp/reference.fa.sa

    cd /output/svaba-tmp
    svaba run -t $BAMT -n $BAMN -p $NCORES -a $ID -G reference.fa -z $STRUCTURAL_SVABAARGS
    ## variant vcf
    mv $ID.svaba.{somatic,germline}.* /output
    ## contig bam
    samtools sort -@ $NCORES -m 8G $ID.contigs.bam -o $ID.svaba.contigs.bam
    samtools index -@ $NCORES $ID.svaba.contigs.bam
    mv $ID.svaba.contigs.bam /output
    ## debug
    if [ "$STRUCTURAL_DEBUG" = true ]; then
        mkdir -p /output/debug
        mv $ID.alignments.txt.gz /output/debug/$ID.svaba.alignments.txt.gz
        mv $ID.bps.txt.gz /output/debug/$ID.svaba.bps.txt.gz
        mv $ID.discordant.txt.gz /output/debug/$ID.svaba.discordant.txt.gz
        mv $ID.svaba.unfiltered.* /output/debug
    fi
    ## clean-up
    rm -rf /output/svaba-tmp
    cd -

fi

## Run GRIDSS
if [ "$STRUCTURAL_GRIDSS" = true ]; then
    echo $(date) "Running GRIDSS Somatic" >> /job/docker_script.log

    ## prepare
    export PATH=$PATH:$SENTIEON_INSTALL_DIR/bin
    mkdir -p /output/gridss-tmp
    ln -sfn $ALIGN_FASTA /output/gridss-tmp/reference.fa
    ln -sfn $ALIGN_FASTA.fai /output/gridss-tmp/reference.fa.fai
    ln -sfn $ALIGN_INDEX.dict /output/gridss-tmp/reference.fa.dict
    ln -sfn $ALIGN_INDEX.alt /output/gridss-tmp/reference.fa.alt
    ln -sfn $ALIGN_INDEX.ann /output/gridss-tmp/reference.fa.ann
    ln -sfn $ALIGN_INDEX.pac /output/gridss-tmp/reference.fa.pac
    ln -sfn $ALIGN_INDEX.amb /output/gridss-tmp/reference.fa.amb
    ln -sfn $ALIGN_INDEX.bwt /output/gridss-tmp/reference.fa.bwt
    ln -sfn $ALIGN_INDEX.sa /output/gridss-tmp/reference.fa.sa

    ## GRIDSS
    gridss --jar /opt/gridss/gridss-jar-with-dependencies.jar \
              --reference /output/gridss-tmp/reference.fa \
              --output $PFX-structural-gridss-temp.vcf.gz \
              --assembly $PFX-structural-gridss.bam \
              --threads $NCORES \
              --workingdir /output/gridss-tmp \
              $STRUCTURAL_GRIDSSARGS \
              --labels $BAMT_SAMPLE,$BAMN_SAMPLE \
              $BAMT $BAMN &>> /job/docker_script.log
    samtools index $PFX-structural-gridss.bam

    ## RepeatMasker
    gridss_annotate_vcf_repeatmasker --jar /opt/gridss/gridss-jar-with-dependencies.jar \
                                     --threads $NCORES \
                                     --workingdir /tmp \
                                     --output $PFX-structural-gridss.vcf \
                                     $PFX-structural-gridss-temp.vcf.gz &>> /job/docker_script.log
    rm $PFX-structural-gridss-temp.vcf.gz*
    bgzip $PFX-structural-gridss.vcf
    tabix $PFX-structural-gridss.vcf.gz

    ## GRIPSS (Somatic)
    java -Xms4G -Xmx16G -jar /opt/gridss/gripss.jar \
         -sample $BAMT_SAMPLE \
         -reference $BAMN_SAMPLE \
         -ref_genome /output/gridss-tmp/reference.fa \
         $STRUCTURAL_GRIPSSARGS \
         -vcf $PFX-structural-gridss.vcf.gz \
         -output_dir /output \
         -output_id $ID &>> /job/docker_script.log
    mv /output/$BAMT_SAMPLE.gripss.filtered.$ID.vcf.gz $PFX-structural-gridss-somatic.vcf.gz
    mv /output/$BAMT_SAMPLE.gripss.filtered.$ID.vcf.gz.tbi $PFX-structural-gridss-somatic.vcf.gz.tbi
    mv /output/$BAMT_SAMPLE.gripss.$ID.vcf.gz $PFX-structural-gridss-unfiltered.vcf.gz
    mv /output/$BAMT_SAMPLE.gripss.$ID.vcf.gz.tbi $PFX-structural-gridss-unfiltered.vcf.gz.tbi
    
    ## clean-up
    rm -rf /output/*idx
    rm -rf /output/gridss.tmp*
    rm -rf /output/gridss-tmp

fi

echo $(date) >> /job/docker_done
