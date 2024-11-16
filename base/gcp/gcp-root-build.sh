#!/usr/bin/env bash
set -eu

GCP_TPO_BOOT_IMAGE=tpoboot-$(echo $TPO_BOOT_VER | tr "." "-")
GCP_TPO_ROOT_IMAGE=tporoot-$(echo $TPO_ROOT_VER | tr "." "-")
GCP_TPO_WORK_IMAGE=tpowork-$(echo $GCP_WORK_VOLUME | tr "." "-")
gcloud compute instances create $GCP_NAME \
    --project=$GCP_PROJECT \
    --zone=$GCP_ZONE \
    --machine-type=n2-standard-32 \
    --subnet=$GCP_SUBNET \
    --network-tier=PREMIUM \
    --tags=$GCP_TAG \
    --maintenance-policy=MIGRATE \
    --scopes "https://www.googleapis.com/auth/cloud-platform" \
    --service-account=$GCP_SERVICE_ACCOUNT \
    --image=$GCP_TPO_BOOT_IMAGE \
    --boot-disk-size=$GCP_BOOT_SIZE --boot-disk-type=pd-standard \
    --create-disk=mode=rw,size=$GCP_ROOT_SIZE,type=pd-ssd,auto-delete=no,device-name=tpo-root,mode=rw,name=$GCP_TPO_ROOT_IMAGE \
    --create-disk=mode=rw,size=1000,type=pd-ssd,auto-delete=yes,device-name=tpo-build,mode=rw,name=$GCP_TPO_ROOT_IMAGE-build \
    --metadata ^@^enable-oslogin=TRUE@startup-script="
#!/usr/bin/env bash
set -eu
sleep 20

## Secrets
export SENTIEON_LICENSE=$SECRETS_SENTIEON_LICENSE

## import common functions
source /usr/bin/bash_funcs.sh

## format disks
mkfs.ext4 -m 0 -F -E lazy_itable_init=0,lazy_journal_init=0,discard /dev/disk/by-id/google-tpo-root
mkdir -p /mnt/disks/tpo-root
mount -o discard,defaults /dev/disk/by-id/google-tpo-root /mnt/disks/tpo-root
chmod a+w /mnt/disks/tpo-root

mkfs.ext4 -m 0 -F -E lazy_itable_init=0,lazy_journal_init=0,discard /dev/disk/by-id/google-tpo-build
mkdir -p /mnt/disks/tpo-build
mount -o discard,defaults /dev/disk/by-id/google-tpo-build /mnt/disks/tpo-build
chmod a+w /mnt/disks/tpo-build

## pull data
mkdir -p /mnt/disks/tpo-build/tpo/{indices,tmp}
chmod -R 777 /tmp
mkdir -p /mnt/disks/tpo-build/refs
gsutil_fast rsync -r gs://$TPO_REFS_VER /mnt/disks/tpo-build/refs
mv /mnt/disks/tpo-build/refs/refs /mnt/disks/tpo-build/tpo/refs

## setup code
mkdir -p /tmp/code/base/root
gsutil_fast rsync -r gs://$RUNTIME_WORK_BUCKET/root /tmp/code/base/root
chmod -R 777 /tmp/code

NCORES=\$(nproc --all)

## BUILD REFS
docker run --privileged --rm=true \
    -e SENTIEON_LICENSE \
    -e NCORES=\$NCORES \
    -v /tmp/code:/code \
    -e CORDS_ALIGN_FASTA=$CORDS_ALIGN_FASTA \
    -e CORDS_ALIGN_NAME=$CORDS_ALIGN_NAME \
    -e CORDS_ALIGN_MASK=$CORDS_ALIGN_MASK \
    -e CRISP_ALIGN_FASTA=$CRISP_ALIGN_FASTA \
    -e CRISP_ALIGN_NAME=$CRISP_ALIGN_NAME \
    -e CRISP_ALIGN_MASK=$CRISP_ALIGN_MASK \
    -e CRISP_QUANT_FASTA=$CRISP_QUANT_FASTA \
    -e CRISP_QUANT_NAME=$CRISP_QUANT_NAME \
    -e CRISP_QUANT_GTF=$CRISP_QUANT_GTF \
    -e CARAT_ANNO_VEPCACHE=$CARAT_ANNO_VEPCACHE \
    -v /mnt/disks/tpo-build/tpo:/tpo \
    -v /mnt/disks/tpo-build/refs:/input \
    gcr.io/$GCP_PROJECT/tporefs:$TPO_BOOT_VER /code/base/root/build_refs.sh &

## push tpo to final disk
wait
sleep 60
rsync -a /mnt/disks/tpo-build/tpo/ /mnt/disks/tpo-root/

## shutdown
sleep 60
retry 3 umount /mnt/disks/tpo-build
retry 3 umount /mnt/disks/tpo-root
retry 3 gcloud compute instances delete -q --project=$GCP_PROJECT --zone=$GCP_ZONE $GCP_NAME

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
gcloud compute images create $GCP_TPO_ROOT_IMAGE --source-disk=$GCP_TPO_ROOT_IMAGE --project=$GCP_PROJECT --source-disk-zone=$GCP_ZONE

## delete disk
sleep 30
gcloud compute disks delete -q $GCP_TPO_ROOT_IMAGE --zone=$GCP_ZONE --project=$GCP_PROJECT
