#!/usr/bin/env bash
set -eu

GCP_TPO_WORK_IMAGE=tpowork-$(echo $GCP_WORK_VOLUME | tr "." "-")
gcloud compute instances create $GCP_NAME \
    --project=$GCP_PROJECT \
    --zone=$GCP_ZONE \
    --machine-type=n1-standard-8 \
    --subnet=$GCP_SUBNET \
    --network-tier=PREMIUM \
    --tags=$GCP_TAG \
    --maintenance-policy=MIGRATE \
    --scopes "https://www.googleapis.com/auth/cloud-platform" \
    --service-account=$GCP_SERVICE_ACCOUNT \
    $GCP_IMAGE \
    --boot-disk-size=50GB --boot-disk-type=pd-standard \
    --create-disk=mode=rw,size=$GCP_WORK_SIZE,type=pd-ssd,auto-delete=no,device-name=tpo-work,mode=rw,name=$GCP_TPO_WORK_IMAGE \
    --metadata ^@^enable-oslogin=TRUE@startup-script="
#!/usr/bin/env bash
set -eu
sleep 20

## format work disk
mkfs.ext4 -m 0 -F -E lazy_itable_init=0,lazy_journal_init=0,discard /dev/disk/by-id/google-tpo-work
mkdir -p /mnt/disks/tpo-work
mount -o discard,defaults /dev/disk/by-id/google-tpo-work /mnt/disks/tpo-work
chmod a+w /mnt/disks/tpo-work

## shutdown
cd /
sleep 60
umount /mnt/disks/tpo-work
gcloud compute instances delete -q --project=$GCP_PROJECT --zone=$GCP_ZONE $GCP_NAME

"

## wait until cloud-build finishes
set +e
while true ; do
    INSTANCE_RUNNIG=$(gcloud compute instances list --project=$GCP_PROJECT --filter="NAME=$GCP_NAME" --format='value(NAME)')
    echo $INSTANCE_RUNNIG
    if [ -z "$INSTANCE_RUNNIG" ]; then
        break
    fi
    sleep 30
done
set -e

## create image from disk
sleep 30
gcloud compute images create $GCP_TPO_WORK_IMAGE --source-disk=$GCP_TPO_WORK_IMAGE --project=$GCP_PROJECT --source-disk-zone=$GCP_ZONE

## delete disk
sleep 30
gcloud compute disks delete -q $GCP_TPO_WORK_IMAGE --zone=$GCP_ZONE --project=$GCP_PROJECT
