#!/usr/bin/env bash
set -eu

## Secrets
export AWS_ACCESS_KEY_ID=$SECRETS_AWS_ACCESS_KEY_ID
export AWS_SECRET_ACCESS_KEY=$SECRETS_AWS_SECRET_ACCESS_KEY
export SENTIEON_LICENSE=$SECRETS_SENTIEON_LICENSE 

## import common functions
source $ROOT/pipe/common/bash_funcs.sh

## setup work directories
JOB=$RUNTIME_RUNS/crisp-align/$GCP_NAME
mkdir -p $JOB/{tmp,input,output}

## setup run folder
source $JOB/config.txt

## import data
retry 3 gsutil_fast cp -n $FQ1 $JOB/input
retry 3 gsutil_fast cp -n $FQ2 $JOB/input

docker run --rm=true \
       -e SENTIEON_LICENSE \
       -v $ROOT:/code \
       -v $RUNTIME_ROOT:/tpo \
       -v $JOB/tmp:/tmp \
       -v $JOB:/job \
       -v $JOB/input:/input \
       -v $JOB/output:/output \
       gcr.io/$GCP_PROJECT/tpocrisp:$TPO_BOOT_VER /code/pipe/crisp-align/crisp-align_docker.sh
