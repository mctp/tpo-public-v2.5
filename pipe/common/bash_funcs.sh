#!/usr/bin/env bash

function gsutil_fast {
    NCORES=$(nproc --all)
    PC=$(expr $NCORES / 2)
    gsutil -m -q \
        -o GSUtil:parallel_process_count=$PC \
        -o GSUtil:parallel_thread_count=2 \
        -o GSUtil:parallel_composite_upload_threshold=100M \
        "$@"
}

function retry {
    n=0
    echo "$2 ${@:3}"
    until [ $n -ge "$1" ]
    do
        $2 "${@:3}" && break
        n=$[$n+1]
        echo "failed... $n"
        sleep 10
    done
}

alias get_startup_script='curl "http://metadata.google.internal/computeMetadata/v1/instance/attributes/startup-script" -H "Metadata-Flavor: Google"'
alias get_startup_output='sudo journalctl -u google-startup-scripts.service '
