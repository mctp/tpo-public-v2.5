#!/bin/bash
# This script initialized a new database with the name provided on an existing GCP instance found using the CURVE_DB variable
set -eu

NEW_CURVE_NAME=$1
TPO_BASE=$2

if [[ -n "$CURVE_DB" ]]; then

  ### initialize the new database on the existing GCP curve instance ###
  args=$(echo $CURVE_DB | sed 's_:_ _g' | xargs printf "host=%s port=%s dbname=%s user=%s password=%s")
  cmd="psql \"$args\""

  #generate the CREATE command
  eval "$cmd -c \"CREATE DATABASE $NEW_CURVE_NAME;\""

  #generate the SCHEMA command and FUNCTIONS
  if [[ -e "$TPO_BASE/rlibs/curve/exec/create_db.sql" ]]; then
    args2=$(echo $CURVE_DB | sed 's_:_ _g' | cut -f1,2,4,5 -d' ' | xargs printf "host=%s port=%s dbname=$NEW_CURVE_NAME user=%s password=%s")
    cmd2="psql \"$args2\""
    eval "$cmd2 < $TPO_BASE/rlibs/curve/exec/create_db.sql"
    eval "$cmd2 < $TPO_BASE/rlibs/curve/exec/functions.sql"
  fi

  # get the name of the existing GCP instance for gcloud commands
  CURVE_IP=$(echo $CURVE_DB | cut -f1 -d:)
  INST_NAME=$(gcloud compute instances list | grep "$CURVE_IP " | cut -f1 -d' ')
  echo "setting up DB on $INST_NAME"

  # generate the command to refresh views in the new database and put it in the crontab
  refresh_cmd="(crontab -l; echo '@daily psql \"host=localhost port=5432 dbname=$NEW_CURVE_NAME user=$SECRETS_CURVE_USER password=$SECRETS_CURVE_PW\" -c \"refresh materialized view somatic_recur;\" -c \"refresh materialized view germline_recur;\" -c \"refresh materialized view fusion_recur;\"')  | crontab -u root -"
  
  # generate the commands to initialize the gene expression table using the GTF file
  pull_cmd="docker pull gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER && docker pull gcr.io/$GCP_PROJECT/tpocode:$TPO_CODE_VER"
  export_cmd="export CURVE_DB=\"\$(hostname -I | cut -f1 -d' '):5432:$NEW_CURVE_NAME:$SECRETS_CURVE_USER:$SECRETS_CURVE_PW\""
  
  ### Make a shell script of the necesary commands, put it on the instance and execute the command via SSH ###
  echo $refresh_cmd > /tmp/curve_db_add_script.sh
  echo $pull_cmd >> /tmp/curve_db_add_script.sh
  echo $export_cmd >> /tmp/curve_db_add_script.sh

  gcloud compute --project $GCP_PROJECT scp --internal-ip /tmp/curve_db_add_script.sh $INST_NAME:/work/
  gcloud compute --project $GCP_PROJECT ssh --internal-ip $INST_NAME --command \
    'chmod o+x /work/curve_db_add_script.sh && sudo $_ && mv /work/curve_db_add_script.sh /work/curve_db_add_script.sh.bak.$(date -I)'

  #clean up 
  rm /tmp/curve_db_add_script.sh
else
  echo "please set CURVE_DB env variable"
fi
