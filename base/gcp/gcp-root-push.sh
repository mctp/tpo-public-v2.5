#!/usr/bin/env bash
set -eu

GCP_TPO_BOOT_IMAGE=tpoboot-$(echo $TPO_BOOT_VER | tr "." "-")
GCP_TPO_ROOT_IMAGE=tporoot-$(echo $TPO_ROOT_VER | tr "." "-")
GCP_TPO_WORK_IMAGE=tpowork-$(echo $GCP_WORK_VOLUME | tr "." "-")
gcloud compute instances create $GCP_NAME \
    --project=$GCP_PROJECT \
    --zone=$GCP_ZONE \
    --machine-type=n2-highcpu-16 \
    --subnet=$GCP_SUBNET \
    --network-tier=PREMIUM \
    --tags=$GCP_TAG \
    --scopes "https://www.googleapis.com/auth/cloud-platform" \
    --service-account=$GCP_SERVICE_ACCOUNT \
    --image=$GCP_TPO_BOOT_IMAGE \
    --boot-disk-size=$GCP_BOOT_SIZE --boot-disk-type=pd-standard --boot-disk-device-name=boot-dev \
    --create-disk=mode=rw,size=$GCP_WORK_DISK,type=pd-ssd,device-name=tpo-work,name=$GCP_NAME-tpo-work,image=$GCP_TPO_WORK_IMAGE  \
    --create-disk=mode=rw,size=$GCP_ROOT_SIZE,type=pd-ssd,device-name=tpo-root,name=$GCP_NAME-tpo-root,image=$GCP_TPO_ROOT_IMAGE \
    --metadata ^@^startup-script="
#!/usr/bin/env bash
set -eu
sleep 20

source /usr/bin/bash_funcs.sh

cd /
tar -cf - tpo/ | pigz -p 32 > /work/tporoot_$TPO_ROOT_VER.tar.gz
retry 3 gsutil_fast cp /work/tporoot_$TPO_ROOT_VER.tar.gz gs://$RUNTIME_WORK_BUCKET/tporoot_$TPO_ROOT_VER.tar.gz
sleep 30
retry 3 gcloud compute instances delete -q --project=$GCP_PROJECT --zone=$GCP_ZONE $GCP_NAME

"@shutdown-script="
#!/usr/bin/env bash
set -eu

"@enable-oslogin=TRUE@user-data='

#cloud-config

bootcmd:
- mkdir -p /tpo
- mount -t ext4 /dev/disk/by-id/google-tpo-root /tpo
- resize2fs /dev/disk/by-id/google-tpo-root
- mkdir -p /work
- mount -t ext4 /dev/disk/by-id/google-tpo-work /work
- resize2fs /dev/disk/by-id/google-tpo-work

'
