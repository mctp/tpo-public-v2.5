#!/usr/bin/env bash
set -eu

cp $RUNTIME_REFS/tools/gridss-2.13.2.tar.gz $ROOT/base/images/tpocords/gridss.tar.gz
cp $RUNTIME_REFS/tools/manta-1.6.0.centos6_x86_64.tar.bz2 $ROOT/base/images/tpocords/manta.tar.bz2
cp $RUNTIME_REFS/tools/strelka-2.9.10.centos6_x86_64.tar.bz2 $ROOT/base/images/tpocords/strelka.tar.bz2
cp $RUNTIME_REFS/tools/biobambam2-2.0.146-release-20191030105216-x86_64-linux-gnu.tar.xz $ROOT/base/images/tpocords/biobambam2-x86_64-linux-gnu.tar.xz
cp $RUNTIME_REFS/tools/gripss-2.0.jar $ROOT/base/images/tpocords/gripss.jar

docker tag gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER tpobase:temp
docker build $BUILD_DOCKER_NOCACHE -t gcr.io/$GCP_PROJECT/tpocords:$TPO_BOOT_VER $ROOT/base/images/tpocords
docker rmi tpobase:temp

rm $ROOT/base/images/tpocords/gridss.tar.gz
rm $ROOT/base/images/tpocords/manta.tar.bz2
rm $ROOT/base/images/tpocords/strelka.tar.bz2
rm $ROOT/base/images/tpocords/biobambam2-x86_64-linux-gnu.tar.xz
rm $ROOT/base/images/tpocords/gripss.jar
