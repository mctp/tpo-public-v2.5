#!/usr/bin/env bash
set -eu 

cp $RUNTIME_REFS/tools/inchworm_19032017.tar.gz $ROOT/base/images/tpocrisp/inchworm.tar.gz
cp $RUNTIME_REFS/tools/mixcr-3.0.13.zip $ROOT/base/images/tpocrisp/mixcr-3.0.13.zip

docker tag gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER tpobase:temp
docker build $BUILD_DOCKER_NOCACHE -t gcr.io/$GCP_PROJECT/tpocrisp:$TPO_BOOT_VER $ROOT/base/images/tpocrisp
docker rmi tpobase:temp

rm $ROOT/base/images/tpocrisp/inchworm.tar.gz
rm $ROOT/base/images/tpocrisp/mixcr-3.0.13.zip
