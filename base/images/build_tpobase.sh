#!/usr/bin/env bash
set -eu

## build docker image
cp $RUNTIME_REFS/tools/gmap-gsnap-2017-11-15.tar.gz $ROOT/base/images/tpobase/gmap-gsnap.tar.gz
cp $RUNTIME_REFS/tools/bcl2fastq2-v2.20.0.422-Linux-x86_64_gcc7.3.deb $ROOT/base/images/tpobase/bcl2fastq2.deb
cp $RUNTIME_REFS/tools/sentieon-genomics-202112.06.tar.gz $ROOT/base/images/tpobase/sentieon-genomics.tar.gz
cp $RUNTIME_REFS/tools/BBMap_37.36.tar.gz $ROOT/base/images/tpobase/bbmap.tar.gz
cp $RUNTIME_REFS/tools/msisensor2-v0.1 $ROOT/base/images/tpobase/msisensor2

docker build $BUILD_DOCKER_NOCACHE \
       -t gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER $ROOT/base/images/tpobase

rm $ROOT/base/images/tpobase/bcl2fastq2.deb
rm $ROOT/base/images/tpobase/sentieon-genomics.tar.gz
rm $ROOT/base/images/tpobase/gmap-gsnap.tar.gz
rm $ROOT/base/images/tpobase/msisensor2
rm $ROOT/base/images/tpobase/bbmap.tar.gz
