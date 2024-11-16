#!/usr/bin/env bash
set -eu

## Secrets
export AWS_ACCESS_KEY_ID=$SECRETS_AWS_ACCESS_KEY_ID
export AWS_SECRET_ACCESS_KEY=$SECRETS_AWS_SECRET_ACCESS_KEY
export SENTIEON_LICENSE=$SECRETS_SENTIEON_LICENSE 

## import common functions
source $ROOT/pipe/common/bash_funcs.sh

## setup work directories
JOB=$RUNTIME_RUNS/bcl/$GCP_NAME
mkdir -p $JOB/{tmp,input,output}
mkdir -p $JOB/input/flowcell

## setup run folder
source $JOB/config.txt

## get the barcode sheet if needed
if [[ $BCL2FASTQ_SHEET =~ gs://* ]]; then
  retry 3 gsutil_fast cp $BCL2FASTQ_SHEET $JOB/input/custom_barcode.csv
fi

TARFILE=$(basename $TAR)
FCID=${TARFILE%.*}

## TAR can be a TARFILE or a directory when running in local mode
if [[ $TAR == gs://* ]] || [[ $TAR == s3://* ]]; then
    retry 3 gsutil_fast cp -n $TAR $JOB/input/
    ## untar sequencing run archive
    tar -xf $JOB/input/$TARFILE -C $JOB/input/flowcell
    ## compatibility with old tar archive structure
    if [ ! -d $JOB/input/flowcell/$FCID ]; then
        mv $JOB/input/flowcell/Illumina/$FCID $JOB/input/flowcell
        rm -rf $JOB/input/flowcell/Illumina
    fi
    FCPATH=$JOB/input/flowcell/$FCID
else
    FCPATH=$TAR
fi

if [[ $LIB == gs://* ]] || [[ $LIB == s3://* ]]; then
    retry 3 gsutil_fast cp -n $LIB $JOB/input/
else
    cp $LIB $JOB/input/
fi

docker run --privileged=true --rm=true \
       -e SENTIEON_LICENSE \
       -v $ROOT:/code \
       -v $RUNTIME_ROOT:/tpo \
       -v $JOB/tmp:/tmp \
       -v $JOB:/job \
       -v $JOB/input:/input \
       -v $JOB/output:/output \
       -v $FCPATH:/flowcell:ro \
       gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER /code/pipe/bcl/bcl_bcl2fastq_docker.sh

docker run --rm=true \
       -e SENTIEON_LICENSE \
       -v $ROOT:/code \
       -v $RUNTIME_ROOT:/tpo \
       -v $JOB/tmp:/tmp \
       -v $JOB:/job \
       -v $JOB/input:/input \
       -v $JOB/output:/output \
       gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER /code/pipe/bcl/bcl_rename_docker.sh
