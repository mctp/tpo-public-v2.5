#!/usr/bin/env bash
set -eu

## Secrets
export AWS_ACCESS_KEY_ID=$SECRETS_AWS_ACCESS_KEY_ID
export AWS_SECRET_ACCESS_KEY=$SECRETS_AWS_SECRET_ACCESS_KEY
export SENTIEON_LICENSE=$SECRETS_SENTIEON_LICENSE 

## import common functions
source $ROOT/pipe/common/bash_funcs.sh

## setup work directories
JOB=$RUNTIME_RUNS/cargo-vault/$GCP_NAME
mkdir -p $JOB/{tmp,input,output}

## setup run folder
source $JOB/config.txt

## import data
mkdir -p $JOB/input/carat-anno
mkdir -p $JOB/input/cords-cnvex
mkdir -p $JOB/input/crisp-tquasr
mkdir -p $JOB/input/crisp-nquasr
mkdir -p $JOB/input/crisp-codac
if [ ! -z "$CARAT_ANNO" ]; then
   gsutil_fast rsync $CARAT_ANNO $JOB/input/carat-anno
fi
if [ ! -z "$CORDS_CNVEX" ]; then
   gsutil_fast rsync $CORDS_CNVEX $JOB/input/cords-cnvex
fi
if [ ! -z "$CRISP_QUASR_TUMOR" ]; then
   gsutil_fast rsync $CRISP_QUASR_TUMOR $JOB/input/crisp-tquasr
fi
if [ ! -z "$CRISP_QUASR_NORMAL" ]; then
   gsutil_fast rsync $CRISP_QUASR_NORMAL $JOB/input/crisp-nquasr
fi
if [ ! -z "$CRISP_CODAC" ]; then
   gsutil_fast rsync $CRISP_CODAC $JOB/input/crisp-codac
fi
if [[ $VAULT_SETTINGS == gs://* ]]; then
    VAULT_SETTINGS_FILE=$(basename $VAULT_SETTINGS)
    retry 3 gsutil_fast cp $VAULT_SETTINGS $JOB/input/$VAULT_SETTINGS_FILE
fi

docker run --rm=true \
       -e SENTIEON_LICENSE \
       -v $ROOT:/code \
       -v $RUNTIME_ROOT:/tpo \
       -v $JOB/tmp:/tmp \
       -v $JOB:/job \
       -v $JOB/input:/input \
       -v $JOB/output:/output \
       gcr.io/$GCP_PROJECT/tpocarat:$TPO_BOOT_VER /code/pipe/cargo-vault/cargo-vault_docker.sh
