#!/usr/bin/env bash
set -eu

docker tag gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER tpobase:temp
docker build $BUILD_DOCKER_NOCACHE -t gcr.io/$GCP_PROJECT/tpocarat:$TPO_BOOT_VER $ROOT/base/images/tpocarat
docker rmi tpobase:temp
