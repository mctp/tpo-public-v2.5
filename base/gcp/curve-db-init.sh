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
    --scopes "https://www.googleapis.com/auth/cloud-platform" \
    --service-account=$GCP_SERVICE_ACCOUNT \
    --image=$GCP_TPO_BOOT_IMAGE \
    --boot-disk-size=$GCP_BOOT_SIZE --boot-disk-type=pd-standard --boot-disk-device-name=boot-dev \
    --create-disk=mode=rw,size=$GCP_WORK_DISK,type=pd-ssd,device-name=tpo-work,name=$GCP_NAME-tpo-work,image=$GCP_TPO_WORK_IMAGE \
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


## import common functions
source /usr/bin/bash_funcs.sh

## install cron
apt -y install cron

# Install postgres and configure
echo \"listen_addresses = '*'\" >> /etc/postgresql/14/main/postgresql.conf
sed -i \"s/peer$/trust/g\" /etc/postgresql/14/main/pg_hba.conf
sed -i \"s/scram-sha-256$/md5/g\" /etc/postgresql/14/main/pg_hba.conf
echo \"host all all 0.0.0.0/0 md5\" >> /etc/postgresql/14/main/pg_hba.conf
echo \"local all all md5\" >> /etc/postgresql/14/main/pg_hba.conf
echo \"password_encryption = md5\" >> /etc/postgresql/14/main/postgresql.conf
sudo service postgresql restart

# Add the user and grant
sudo -u postgres psql << EOF 
CREATE USER $SECRETS_CURVE_USER WITH PASSWORD '$SECRETS_CURVE_PW';
CREATE database curvedb;
GRANT CONNECT ON DATABASE curvedb TO $SECRETS_CURVE_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public to $SECRETS_CURVE_USER;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA public to $SECRETS_CURVE_USER;
ALTER USER $SECRETS_CURVE_USER CREATEDB;
EOF
sed -i \"s/trust$/peer/g\" /etc/postgresql/14/main/pg_hba.conf

## extract the schema from the code docker image
gcloud auth configure-docker --quiet
docker pull gcr.io/$GCP_PROJECT/tpocode:$TPO_CODE_VER
CODE=\$(docker run -d \
       -v /work:/work \
       gcr.io/$GCP_PROJECT/tpocode:$TPO_CODE_VER cp -r /code/rlibs/curve/exec /work)

# create the schema and functions
psql \"host=localhost port=5432 dbname=curvedb user=$SECRETS_CURVE_USER password=$SECRETS_CURVE_PW\" < /work/exec/create_db.sql
psql \"host=localhost port=5432 dbname=curvedb user=$SECRETS_CURVE_USER password=$SECRETS_CURVE_PW\" < /work/exec/functions.sql

# set up the refresh
echo '0 0 * * * psql \"host=localhost port=5432 dbname=curvedb user=$SECRETS_CURVE_USER password=$SECRETS_CURVE_PW\" -c \"refresh materialized view somatic_recur;\" -c \"refresh materialized view germline_recur;\" -c \"refresh materialized view fusion_recur;\"'  | crontab -



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
