#!/usr/bin/env bash
set -eu

export HOME=~
gcloud auth configure-docker --quiet
docker pull gcr.io/$GCP_PROJECT/tpocode:$TPO_CODE_VER &>> /dev/null
