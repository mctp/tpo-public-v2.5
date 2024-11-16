#!/usr/bin/env bash
set -eu

## create the schema and functions
DB_CONN="host=localhost port=$CURVE_DB_PORT dbname=$CURVE_DB_NAME user=$CURVE_DB_USER password=$CURVE_DB_PASS"
sudo -u postgres psql -c "CREATE DATABASE $CURVE_DB_NAME;"
sudo psql "$DB_CONN" < $ROOT/create_db.sql
sudo psql "$DB_CONN" < $ROOT/functions.sql

## set up the refresh
echo "@daily psql \"$DB_CONN\" -c \"refresh materialized view somatic_recur;\" -c \"refresh materialized view germline_recur;\" -c \"refresh materialized view fusion_recur;\"" | sudo crontab -u root -
