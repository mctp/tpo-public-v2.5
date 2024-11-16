#!/usr/bin/env bash
set -eu

docker tag gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER tpobase:temp
docker build $BUILD_DOCKER_NOCACHE -t gcr.io/$GCP_PROJECT/tporefs:$TPO_BOOT_VER $ROOT/base/images/tporefs
docker rmi tpobase:temp
