#!/bin/bash
set -eu

## copy scripts and schema
gcloud compute --project $GCP_PROJECT scp --zone $GCP_ZONE --internal-ip $ROOT/rlibs/curve/exec/curve-add-db.sh $CURVE_INSTANCE:/tmp
gcloud compute --project $GCP_PROJECT scp --zone $GCP_ZONE --internal-ip $ROOT/rlibs/curve/exec/create_db.sql $CURVE_INSTANCE:/tmp
gcloud compute --project $GCP_PROJECT scp --zone $GCP_ZONE --internal-ip $ROOT/rlibs/curve/exec/functions.sql $CURVE_INSTANCE:/tmp

## add database in host
gcloud compute --project $GCP_PROJECT ssh --zone $GCP_ZONE --internal-ip $CURVE_INSTANCE --command \
  "export ROOT=/tmp CURVE_DB_PORT=$CURVE_DB_PORT CURVE_DB_NAME=$CURVE_DB_NAME CURVE_DB_USER=$CURVE_DB_USER CURVE_DB_PASS=$CURVE_DB_PASS && bash /tmp/curve-add-db.sh"
