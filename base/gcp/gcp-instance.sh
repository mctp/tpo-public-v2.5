#!/usr/bin/env bash
set -eu

GCP_TPO_BOOT_IMAGE=tpoboot-$(echo $TPO_BOOT_VER | tr "." "-")
GCP_TPO_ROOT_IMAGE=tporoot-$(echo $TPO_ROOT_VER | tr "." "-")
GCP_TPO_WORK_IMAGE=tpowork-$(echo $GCP_WORK_VOLUME | tr "." "-")

if [[ ! -z "$GCP_NFS_SHARE" ]]; then
    GCP_NFS_SHARE_IP=$(gcloud filestore instances describe $GCP_NFS_SHARE \
                              --project=$GCP_PROJECT \
                              --zone=$GCP_ZONE \
                              --format="csv(fileShares[0].name[],networks[0].ipAddresses[0])" \
                           | grep share | cut -f2 -d,)
else
    GCP_NFS_SHARE_IP=""
fi

if [ "$GCP_WORK_DISK_TYPE" == "local-ssd" ]; then
    WORK_DISK_LINE="--local-ssd interface=NVME"
elif [ "$GCP_WORK_DISK_TYPE" == "boot-disk" ]; then
    WORK_DISK_LINE=""
else
    WORK_DISK_LINE="--create-disk=mode=rw,size=$GCP_WORK_DISK,type=$GCP_WORK_DISK_TYPE,device-name=tpo-work,name=$GCP_NAME-tpo-work,image=$GCP_TPO_WORK_IMAGE"
fi

gcloud compute instances create $GCP_NAME \
    --project=$GCP_PROJECT \
    --zone=$GCP_ZONE \
    --machine-type=$GCP_MACHINE \
    --subnet=$GCP_SUBNET \
    --network-tier=PREMIUM \
    --tags=$GCP_TAG \
    --scopes "https://www.googleapis.com/auth/cloud-platform" \
    --service-account=$GCP_SERVICE_ACCOUNT \
    --image=$GCP_TPO_BOOT_IMAGE \
    --boot-disk-size=$GCP_BOOT_DISK --boot-disk-type=pd-standard --boot-disk-device-name=boot-dev \
    $WORK_DISK_LINE \
    --create-disk=mode=rw,size=$GCP_ROOT_SIZE,type=pd-ssd,device-name=tpo-root,name=$GCP_NAME-tpo-root,image=$GCP_TPO_ROOT_IMAGE \
    $GCP_MACHINE_ARGS \
    --metadata ^@^startup-script="
#!/usr/bin/env bash
set -eu
sleep 20

## Secrets
export AWS_ACCESS_KEY_ID=$SECRETS_AWS_ACCESS_KEY_ID
export AWS_SECRET_ACCESS_KEY=$SECRETS_AWS_SECRET_ACCESS_KEY
export SENTIEON_LICENSE=$SECRETS_SENTIEON_LICENSE

if [[ ! -z \"$GCP_NFS_SHARE_IP\" ]]; then
   mkdir -p /mnt/share
   mount $GCP_NFS_SHARE_IP:/share /mnt/share
   chmod 777 /mnt/share
fi

## import common functions
source /usr/bin/bash_funcs.sh

if [[ ! -z \"$GCP_JOB\" ]]; then
   
   ## get script from gs
   mkdir -p /work/job /work/tmp
   retry 3 gsutil_fast rsync -r $GCP_JOB /work/job
   chmod +x /work/job/script.sh
   
   ## run gcp
   cd /work/job
   ./script.sh
   
   ## sync results
   retry 3 gsutil_fast rsync -r /work/job $GCP_JOB

   ## shutdown
   retry 3 gcloud compute instances delete -q --project=$GCP_PROJECT --zone=$GCP_ZONE $GCP_NAME

fi

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
- if [ -e /dev/disk/by-id/google-tpo-work ]; then
      mount -t ext4 /dev/disk/by-id/google-tpo-work /work;
      resize2fs /dev/disk/by-id/google-tpo-work;
  elif [ -e /dev/nvme0n1 ]; then
      mkfs.ext4 /dev/nvme0n1;
      mount /dev/nvme0n1 /work;
  fi

'
