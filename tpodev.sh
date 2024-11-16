#!/usr/bin/env bash

help_msg() {
  echo "Command: tpodev.sh
Usage: 
  tpodev.sh tpo_config container_name [additional mountpoints]"
  exit 1
}

if [[ ! -n "$1" || ! -n "$2" ]]; then
    help_msg
fi

TPO_BOOT_VER=$(sed -nr "/^\[TPO\]/ { :l /^BOOT_VER[ ]*=/ { s/[^=]*=[ ]*//; p; q;}; n; b l;}" $1)
RUNTIME_ROOT=$(sed -nr "/^\[RUNTIME\]/ { :l /^ROOT[ ]*=/ { s/[^=]*=[ ]*//; p; q;}; n; b l;}" $1)
RUNTIME_WORK=$(sed -nr "/^\[RUNTIME\]/ { :l /^WORK[ ]*=/ { s/[^=]*=[ ]*//; p; q;}; n; b l;}" $1)
RUNTIME_TEMP=$(sed -nr "/^\[RUNTIME\]/ { :l /^TEMP[ ]*=/ { s/[^=]*=[ ]*//; p; q;}; n; b l;}" $1)
TPO_CODE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

WHAT_DOCKER=$(docker --version)
case $WHAT_DOCKER in
    Docker*) USER="--user=$(id -u) -w /home/$(id -nu)" ;;
    podman*) USER="" ;;
    * ) echo "No Docker or Podman is found."; exit 0 ;;
esac

SENTIEON_LICENSE=$(sed -nr "/^\[SECRETS\]/ { :l /^SENTIEON_LICENSE[ ]*=/ { s/[^=]*=[ ]*//; p; q;}; n; b l;}" $1)

docker run -t -i --rm \
    "${@:3}" \
    $USER \
    -e SENTIEON_LICENSE=$SENTIEON_LICENSE \
    --name $2 \
    -v $RUNTIME_ROOT:/tpo \
    -v $RUNTIME_WORK:/work \
    -v $RUNTIME_TEMP:/tmp \
    -v $TPO_CODE:/code \
    -v $HOME:/host_home \
    "$(id -nu)-tpodev:$TPO_BOOT_VER" \
    /bin/bash
