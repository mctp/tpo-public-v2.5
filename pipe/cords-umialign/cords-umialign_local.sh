#!/usr/bin/env bash
set -eu

## Secrets
export AWS_ACCESS_KEY_ID=$SECRETS_AWS_ACCESS_KEY_ID
export AWS_SECRET_ACCESS_KEY=$SECRETS_AWS_SECRET_ACCESS_KEY
export SENTIEON_LICENSE=$SECRETS_SENTIEON_LICENSE 

## import common functions
source $ROOT/pipe/common/bash_funcs.sh

## setup work directories
JOB=$RUNTIME_RUNS/cords-umialign/$GCP_NAME
mkdir -p $JOB/{tmp,input,output}

## setup run folder
source $JOB/config.txt

## import data
IFS=';' read -ra FQ1SA <<< $FQ1S
for i in "${FQ1SA[@]}"; do
    retry 3 gsutil_fast cp $i $JOB/input
done
IFS=';' read -ra FQ2SA <<< $FQ2S
for i in "${FQ2SA[@]}"; do
    retry 3 gsutil_fast cp $i $JOB/input
done

docker run --rm=true \
       -e SENTIEON_LICENSE \
       -v $ROOT:/code \
       -v $RUNTIME_ROOT:/tpo \
       -v $JOB/tmp:/tmp \
       -v $JOB:/job \
       -v $JOB/input:/input \
       -v $JOB/output:/output \
       gcr.io/$GCP_PROJECT/tpocords:$TPO_BOOT_VER /code/pipe/cords-umialign/cords-umialign_docker.sh
