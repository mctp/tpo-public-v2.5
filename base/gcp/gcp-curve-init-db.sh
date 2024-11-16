#!/bin/bash
set -eu

## copy curve init script
gcloud compute --project $GCP_PROJECT scp --zone $GCP_ZONE --internal-ip $ROOT/rlibs/curve/exec/curve-init-db.sh $CURVE_INSTANCE:/tmp

## initialize database remotely
gcloud compute --project $GCP_PROJECT ssh --zone $GCP_ZONE --internal-ip $CURVE_INSTANCE --command \
  "export CURVE_DB_USER=$CURVE_DB_USER CURVE_DB_PASS=$CURVE_DB_PASS && bash /tmp/curve-init-db.sh"
