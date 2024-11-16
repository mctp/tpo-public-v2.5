#!/usr/bin/env bash
set -eu

## BUILD CODE
export HOME=~
cd $ROOT


if [ "$BUILD_COMMIT" = "commit" ]; then
    git archive --format=tar.gz -o \
        $ROOT/base/images/tpocode/tpocode.tar.gz HEAD
fi

if [ "$BUILD_COMMIT" = "dirty" ]; then
    uploadStash=`git stash create`; git archive --format=tar.gz -o \
        $ROOT/base/images/tpocode/tpocode.tar.gz ${uploadStash:-HEAD}
fi

docker run --rm=true \
       -v $ROOT:/code \
       gcr.io/$GCP_PROJECT/tpobase:$TPO_BOOT_VER \
       /code/rlibs/install.sh

tar -C $ROOT -cf $ROOT/base/images/tpocode/tporlibs.tar.gz rlibs/install

docker build $BUILD_DOCKER_NOCACHE -t gcr.io/$GCP_PROJECT/tpocode:$TPO_CODE_VER \
       $ROOT/base/images/tpocode

if [ "$BUILD_REMOVE" = "remove" ]; then
    rm -rf $ROOT/rlibs/install/*
fi

rm $ROOT/base/images/tpocode/tpocode.tar.gz
rm $ROOT/base/images/tpocode/tporlibs.tar.gz
