#!/usr/bin/env bash
set -eu
source $ROOT/pipe/common/bash_funcs.sh

if [ -d "$RUNTIME_ROOT" ]; then
    echo "$RUNTIME_ROOT" "exists"
    exit 1
fi

mkdir -p $(dirname $RUNTIME_ROOT)

gsutil_fast cp gs://$TEMP_BUCKET/tporoot_$TPO_ROOT_VER.tar.gz $RUNTIME_TEMP/tporoot_$TPO_ROOT_VER.tar.gz 
tar xf $RUNTIME_TEMP/tporoot_$TPO_ROOT_VER.tar.gz -C $(dirname $RUNTIME_ROOT)
