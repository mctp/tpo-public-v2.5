#!/usr/bin/env bash
set -eu

GCP_TPO_BOOT_IMAGE=tpoboot-$(echo $TPO_BOOT_VER | tr "." "-")
GCP_TPO_ROOT_IMAGE=tporoot-$(echo $TPO_ROOT_VER | tr "." "-")
GCP_TPO_WORK_IMAGE=tpowork-$(echo $GCP_WORK_VOLUME | tr "." "-")
gcloud compute instances create $GCP_NAME \
    --project=$GCP_PROJECT \
    --zone=$GCP_ZONE \
    --machine-type=$GCP_MACHINE \
    --subnet=$GCP_SUBNET \
    --network-tier=PREMIUM \
    --tags=$GCP_TAG \
    --labels=$GCP_LABEL \
    --scopes "https://www.googleapis.com/auth/cloud-platform" \
    --service-account=$GCP_SERVICE_ACCOUNT \
    --image=$GCP_TPO_BOOT_IMAGE \
    --boot-disk-size=$GCP_BOOT_SIZE --boot-disk-type=pd-standard --boot-disk-device-name=boot-dev \
    --create-disk=mode=rw,size=$GCP_WORK_DISK,type=$GCP_WORK_DISK_TYPE,device-name=tpo-work,name=$GCP_NAME-tpo-work,image=$GCP_TPO_WORK_IMAGE \
    --create-disk=mode=rw,size=$GCP_ROOT_SIZE,type=pd-standard,device-name=tpo-root,name=$GCP_NAME-tpo-root,image=$GCP_TPO_ROOT_IMAGE \
    $GCP_MACHINE_ARGS \
    --metadata ^@^startup-script="
#!/usr/bin/env bash
set -eu
sleep 20

## Secrets
export AWS_ACCESS_KEY_ID=$SECRETS_AWS_ACCESS_KEY_ID
export AWS_SECRET_ACCESS_KEY=$SECRETS_AWS_SECRET_ACCESS_KEY
export SENTIEON_LICENSE=$SECRETS_SENTIEON_LICENSE 

## import common functions
source /usr/bin/bash_funcs.sh

## setup work directories
mkdir -p /work/job/{tmp,input,output}
mkdir -p /work/job/input/flowcell
echo \$(date) 'Instance ready' >> /work/job/gcp_script.log

## setup run folder
retry 3 gsutil_fast rsync -r gs://$RUNTIME_WORK_BUCKET/pipe/bcl/$GCP_NAME /work/job
source /work/job/config.txt

TARFILE=\$(basename \$TAR)
FCID=\${TARFILE%.*}

## get the barcode sheet if needed
if [[ \$BCL2FASTQ_SHEET =~ gs://* ]]; then
  retry 3 gsutil_fast cp \$BCL2FASTQ_SHEET /work/job/input/custom_barcode.csv
  cp /work/job/input/custom_barcode.csv /work/job/output/
fi

## import data
retry 3 gsutil_fast cp \$TAR /work/job/input/
retry 3 gsutil_fast cp \$LIB /work/job/input/

## extract data
tar -xf /work/job/input/\$TARFILE -C /work/job/input/flowcell
## compatibility with old tar archive structure
if [ ! -d /work/job/input/flowcell/\$FCID ]; then
    mv /work/job/input/flowcell/Illumina/\$FCID /work/job/input/flowcell
    rm -rf /work/job/input/flowcell/Illumina
fi
FCPATH=/work/job/input/flowcell/\$FCID

## setup test directory
mkdir -p /test
if [ ! -z $TEST_TEST_DRIVE ]; then
   retry 3 gsutil_fast rsync gs://$TEST_TEST_DRIVE /test
fi

## setup code
gcloud auth configure-docker --quiet
docker pull gcr.io/$GCP_PROJECT/tpocode:$TPO_CODE_VER
CODE=\$(docker run -d gcr.io/$GCP_PROJECT/tpocode:$TPO_CODE_VER)

## provenance
cp /work/job/config.txt /work/job/output/\$ID-config.txt
docker images --no-trunc | grep -v '<none>' > /work/job/output/\$ID-docker.txt

## Create FASTQ
docker run --rm=true \
       -e SENTIEON_LICENSE \
       --volumes-from \$CODE \
       -v /tpo:/tpo \
       -v /work/job/tmp:/tmp \
       -v /work/job:/job \
       -v /work/job/input:/input \
       -v /work/job/output:/output \
       -v \$FCPATH:/flowcell:ro \
       -v /test:/test \
       gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER /code/pipe/bcl/bcl_bcl2fastq_docker.sh

## Rename FASTQ
docker run --rm=true \
       -e SENTIEON_LICENSE \
       --volumes-from \$CODE \
       -v /tpo:/tpo \
       -v /work/job/tmp:/tmp \
       -v /work/job:/job \
       -v /work/job/input:/input \
       -v /work/job/output:/output \
       -v /test:/test \
       gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER /code/pipe/bcl/bcl_rename_docker.sh

## wait for docker to finish
while [ ! -f /work/job/docker_done ]; do sleep 10; done
echo \$(date) 'Run finished' >> /work/job/gcp_script.log

retry 3 gsutil_fast rsync -r /work/job/output gs://$RUNTIME_WORK_BUCKET/repo/bcl/\$ID
echo \$(date) 'Run output synced' >> /work/job/gcp_script.log

journalctl -u google-startup-scripts.service --since \"\$(uptime -s)\" > /work/job/gcp_journal.txt

retry 3 gsutil_fast rsync -x 'input|output|docker_start|docker_done' -r /work/job gs://$RUNTIME_WORK_BUCKET/pipe/bcl/$GCP_NAME
echo \$(date) 'Run synced' >> /work/job/gcp_script.log

## shutdown
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
