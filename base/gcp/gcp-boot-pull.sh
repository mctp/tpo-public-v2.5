#!/usr/bin/env bash
set -eu

export HOME=~
gcloud auth configure-docker --quiet
docker pull gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER &>> /dev/null
docker pull gcr.io/$GCP_PROJECT/tporefs:$TPO_BOOT_VER &>> /dev/null
docker pull gcr.io/$GCP_PROJECT/tpocords:$TPO_BOOT_VER &>> /dev/null
docker pull gcr.io/$GCP_PROJECT/tpocrisp:$TPO_BOOT_VER &>> /dev/null
docker pull gcr.io/$GCP_PROJECT/tpocarat:$TPO_BOOT_VER &>> /dev/null
