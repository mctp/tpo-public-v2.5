#!/usr/bin/env bash
set -eu 

docker tag gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER tpobase:temp
docker build $BUILD_DOCKER_NOCACHE -t $(id -nu)-tpodev:$TPO_BOOT_VER \
    --build-arg user=$(id -nu) \
    --build-arg userid=$(id -u) \
    --build-arg passwd=pa55w0rd! \
    $ROOT/base/images/tpodev
docker rmi tpobase:temp
