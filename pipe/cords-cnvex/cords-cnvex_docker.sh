#!/usr/bin/env bash
echo $(date) >> /job/docker_start

set -eu
shopt -s nullglob

#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate cnvex

source /job/config.txt

NCORES=$(nproc --all)


PFX=/output/$ID

## collect input CRAMs/BAMs
echo $(date) "Collecting input BAMs" >> /job/docker_script.log
if [ -z "$ALNT" ]; then
    BAMT=""
else
    BID=$(basename $ALNT)
    if [ -f "/input/$BID/$BID.cram" ]; then
        BAMT="/input/$BID/$BID.cram"
    elif [ -f "/input/$BID/$BID.bam" ]; then
        BAMT="/input/$BID/$BID.bam"
    fi
    BAMT="--tumor $BAMT"
fi
if [ -z "$ALNN" ]; then
    BAMN=""
else
    BID=$(basename $ALNN)
    if [ -f "/input/$BID/$BID.cram" ]; then
        BAMN="/input/$BID/$BID.cram"
    elif [ -f "/input/$BID/$BID.bam" ]; then
        BAMN="/input/$BID/$BID.bam"
    fi
    BAMN="--normal $BAMN"
fi

echo $(date) "Collecting input VCFs" >> /job/docker_script.log
if [ -z "$TVAR" ]; then
    TVCF=""
else
    TVARID=$(basename $TVAR)
    TVCF="--tvcf /input/target/$TVARID/*-dnascope*.vcf.gz"
fi
if [ -z "$GVAR" ]; then
    GVCF=""
else
    GVARID=$(basename $GVAR)
    GVCF="--gvcf /input/genome/$GVARID/*-dnascope*.vcf.gz"
fi
if [ -z "$CNVEX_CAPTURE" ]; then
    CAPT=""
else
    CAPT="--capture $CNVEX_CAPTURE"
fi

echo $(date) "Collecting Pool and Settings" >> /job/docker_script.log
if [ -z "$CNVEX_POOL" ]; then
    POOL=""
else
    POOL="--pool $CNVEX_POOL"
fi
if [ -z "$CNVEX_POPAF" ]; then
    POPAF=""
else
    POPAF="--popaf $CNVEX_POPAF"
fi
if [[ $CNVEX_SETTINGS == gs://* ]]; then
    CNVEX_SETTINGS=/input/$(basename $CNVEX_SETTINGS)
fi
if [ -z "$CNVSEX" ]; then
    SEX=""
else
    SEX="--sex $CNVSEX"
fi

## collect input CNVEX files
if [ -z "$CNVX" ]; then
    CNVX=""
else
    CNVXID=$(basename $CNVX)
    CNVX="--inp /input/$CNVXID/$CNVXID.rds"
    cp /input/$CNVXID/$CNVXID.rds $PFX.rds 2>/dev/null || :
    cp /input/$CNVXID/$CNVXID-opts.rds $PFX-opts.rds 2>/dev/null || :
    cp /input/$CNVXID/$CNVXID-somatic-segment.rds $PFX-somatic-segment.rds 2>/dev/null || :
    cp /input/$CNVXID/$CNVXID-somatic-model.rds $PFX-somatic-model.rds 2>/dev/null || :
    cp /input/$CNVXID/$CNVXID-germline-segment.rds $PFX-germline-segment.rds 2>/dev/null || :
    cp /input/$CNVXID/$CNVXID-germline-model.rds $PFX-germline-model.rds 2>/dev/null || :
fi

if [ "$CNVEX_PROCESS" = true ]; then
    echo $(date) "Starting CNVEX process" >> /job/docker_script.log
    Rscript /code/pipe/cords-cnvex/cnvex_process.R \
        -g $CNVEX_ASSEMBLY \
        -s $CNVEX_SETTINGS \
        -f $ALIGN_FASTA \
        -j 6 \
        $CNVX \
        $CAPT \
        $POOL \
        $POPAF \
        $BAMT \
        $BAMN \
        $TVCF \
        $GVCF \
        $SEX \
        -o $PFX.rds
    echo $(date) "Finished CNVEX process" >> /job/docker_script.log
fi

if [ "$CNVEX_SOMATIC" = true ]; then
    echo $(date) "Starting CNVEX somatic segment" >> /job/docker_script.log
    Rscript /code/pipe/cords-cnvex/cnvex_segment.R \
            -s $CNVEX_SETTINGS \
            -j $NCORES \
            -i $PFX.rds \
            $POOL \
            -a tumor \
            -o $PFX-somatic-segment.rds
    echo $(date) "Finished CNVEX somatic segment" >> /job/docker_script.log
    echo $(date) "Starting CNVEX somatic model search" >> /job/docker_script.log
    Rscript /code/pipe/cords-cnvex/cnvex_search.R \
            -s $CNVEX_SETTINGS \
            -j $NCORES \
            -i $PFX.rds \
            $POOL \
            -a tumor \
            $CNVEX_SOMATICSEARCH \
            -e $PFX-somatic-segment.rds \
            -o $PFX-somatic-model.rds
    echo $(date) "Finished CNVEX somatic model search" >> /job/docker_script.log
fi

if [ "$CNVEX_GERMLINE" = true ]; then
    echo $(date) "Starting CNVEX germline segment" >> /job/docker_script.log
    Rscript /code/pipe/cords-cnvex/cnvex_segment.R \
            -s $CNVEX_SETTINGS \
            -j $NCORES \
            -i $PFX.rds \
            $POOL \
            -a normal \
            -o $PFX-germline-segment.rds
    echo $(date) "Finished CNVEX germline segment" >> /job/docker_script.log
    echo $(date) "Starting CNVEX germline model search" >> /job/docker_script.log
    Rscript /code/pipe/cords-cnvex/cnvex_search.R \
            -s $CNVEX_SETTINGS \
            -j $NCORES \
            -i $PFX.rds \
            $POOL \
            $CNVEX_GERMLINESEARCH \
            -a normal \
            -e $PFX-germline-segment.rds \
            -o $PFX-germline-model.rds
    echo $(date) "Finished CNVEX germline model search" >> /job/docker_script.log
fi

if [ -n "$CNVEX_SOMATICPICK" -a -f "$PFX-somatic-segment.rds" -a -f "$PFX-somatic-model.rds" ]; then
    Rscript /code/pipe/cords-cnvex/cnvex_digest.R \
        -s $CNVEX_SETTINGS \
        -j $NCORES \
        -i $PFX.rds \
        -a tumor \
        $POOL \
        -e $PFX-somatic-segment.rds \
        -m $PFX-somatic-model.rds \
        -p $CNVEX_SOMATICPICK \
        -o $PFX-somatic
fi

if [ -n "$CNVEX_GERMLINEPICK" -a -f "$PFX-germline-segment.rds" -a -f "$PFX-germline-model.rds" ]; then
    Rscript /code/pipe/cords-cnvex/cnvex_digest.R \
        -s $CNVEX_SETTINGS \
        -j $NCORES \
        -i $PFX.rds \
        -a normal \
        $POOL \
        -e $PFX-germline-segment.rds \
        -m $PFX-germline-model.rds \
        -p $CNVEX_GERMLINEPICK \
        -o $PFX-germline
fi

if [ -n "$CNVEX_SOMATICPLOT" ]; then
    if compgen -G "$PFX-somatic-digest*" > /dev/null; then
        ## TODO: support ggopts and override
        Rscript /code/pipe/cords-cnvex/cnvex_ggplots.R \
            -j $NCORES \
            -i $PFX-somatic-digest \
            -o $PFX-somatic \
            -t $CNVEX_SOMATICPLOT
    fi
fi

if [ -n "$CNVEX_GERMLINEPLOT" ]; then
    if compgen -G "$PFX-germline-digest*" > /dev/null; then
        ## TODO: support ggopts and override
        Rscript /code/pipe/cords-cnvex/cnvex_ggplots.R \
            -j $NCORES \
            -i $PFX-germline-digest \
            -o $PFX-germline \
            -t $CNVEX_GERMLINEPLOT
    fi
fi


echo $(date) >> /job/docker_done
