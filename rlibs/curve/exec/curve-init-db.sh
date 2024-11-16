#!/usr/bin/env bash
set -eu

#### Configure PostgreSQL

## Install postgres and configure
echo "listen_addresses = '*'" | sudo tee -a /etc/postgresql/14/main/postgresql.conf
sudo sed -i "s/peer$/trust/g" /etc/postgresql/14/main/pg_hba.conf
sudo sed -i "s/scram-sha-256$/md5/g" /etc/postgresql/14/main/pg_hba.conf
echo "host all all 0.0.0.0/0 md5" | sudo tee -a /etc/postgresql/14/main/pg_hba.conf
echo "local all all md5" | sudo tee -a /etc/postgresql/14/main/pg_hba.conf
echo "password_encryption = md5" | sudo tee -a /etc/postgresql/14/main/postgresql.conf
sudo service postgresql restart

## Add the user and grant permissions
sudo -u postgres psql << EOF
CREATE USER $CURVE_DB_USER WITH PASSWORD '$CURVE_DB_PASS';
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public to $CURVE_DB_USER;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA public to $CURVE_DB_USER;
ALTER USER $CURVE_DB_USER CREATEDB;
EOF
sudo sed -i "s/trust$/peer/g" /etc/postgresql/14/main/pg_hba.conf
sudo service postgresql restart
