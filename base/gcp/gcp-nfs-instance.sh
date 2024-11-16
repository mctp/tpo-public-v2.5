#!/usr/bin/env bash
set -eu

gcloud filestore instances create $GCP_NFS_SHARE \
    --project=$GCP_PROJECT \
    --zone=$GCP_ZONE \
    --tier=$GCP_NFS_TIER \
    --file-share=name=share,capacity=$GCP_NFS_SIZE \
    --network=name=$GCP_NETWORK
