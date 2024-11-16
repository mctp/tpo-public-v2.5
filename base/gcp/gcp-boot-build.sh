#!/usr/bin/env bash
set -eu

GCP_TPO_BOOT_IMAGE=tpoboot-$(echo $TPO_BOOT_VER | tr "." "-")
gcloud compute instances create $GCP_NAME \
    --project=$GCP_PROJECT \
    --zone=$GCP_ZONE \
    --machine-type=n2-standard-8 \
    --subnet=$GCP_SUBNET \
    --network-tier=PREMIUM \
    --tags=$GCP_TAG \
    --maintenance-policy=MIGRATE \
    --scopes "https://www.googleapis.com/auth/cloud-platform" \
    --service-account=$GCP_SERVICE_ACCOUNT \
    $GCP_IMAGE \
    --boot-disk-size=$GCP_BOOT_SIZE --boot-disk-type=pd-standard --boot-disk-device-name=boot-dev --no-boot-disk-auto-delete \
    --metadata ^@^enable-oslogin=TRUE@startup-script="
#!/usr/bin/env bash
set -eu
sleep 60

## add bash helper founctions
gsutil cp gs://$TPO_REFS_VER/tools/bash_funcs.sh /usr/bin/bash_funcs.sh
chmod +x /usr/bin/bash_funcs.sh
source /usr/bin/bash_funcs.sh

## configuration
export DEBIAN_FRONTEND=noninteractive

## update install cleanup
echo 'Update: started'
apt-get update -y
apt-get upgrade -y
echo 'Update: finished update/upgrade'
apt-get install -y apt-transport-https ca-certificates gnupg nfs-common locales tmux software-properties-common docker.io
apt-get install -y python3-pip
apt-get install -y emacs-nox vim-nox pigz ack tree rsync
apt-get install -y postgresql cron
echo 'Update: finished installs'
snap remove google-cloud-sdk
apt remove -y unattended-upgrades
apt autoremove -y
echo 'Update: finished removals'

## enables samtools to access gs://
apt-get install -y libcurl4-gnutls-dev
cd /tmp
wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2
tar xf samtools-1.17.tar.bz2
cd samtools-1.17
./configure --disable-bz2 --disable-lzma --enable-libcurl --without-curses
make -j8
make install
cd -

## normal gcloud/gsutil
echo 'deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main' | sudo tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -
apt-get update -y && apt-get install -y google-cloud-sdk

## locales
locale-gen en_US.UTF-8

## moke
pip3 install moke

## pull docker
gcloud auth configure-docker --quiet
docker pull gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER
docker pull gcr.io/$GCP_PROJECT/tpocords:$TPO_BOOT_VER
docker pull gcr.io/$GCP_PROJECT/tpocrisp:$TPO_BOOT_VER
docker pull gcr.io/$GCP_PROJECT/tpocarat:$TPO_BOOT_VER

sleep 20
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
gcloud compute images create $GCP_TPO_BOOT_IMAGE --source-disk=$GCP_NAME --project=$GCP_PROJECT --source-disk-zone=$GCP_ZONE

## delete disk
gcloud compute disks delete -q $GCP_NAME --zone=$GCP_ZONE --project=$GCP_PROJECT
