#!/usr/bin/env bash
# ShixiangWang@2020
# w_shixiang@163.com
# modified by Ryan Rebernick 2023-03-31

# static path (do NOT change)
GISTIC_LOC=/opt/GISTIC


# Additional info
################
## More about these options, please read 
## ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTICDocumentation_standalone.htm
## https://www.genepattern.org/modules/docs/GISTIC_2.0
## Available reference list
    # hg16.mat/hg17.mat/hg18.mat/hg19.mat
    # hg19.UCSC.add_miR.140312.refgene.mat
    # hg38.UCSC.add_miR.160920.refgene.mat
################


# set params
################
INPUT_FILE=~CUSTOM~
OUTPUT_DIR=~CUSTOM~
REF_PATH=$GISTIC_LOC/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat
################


# Paths
################
mkdir -p $OUTPUT_DIR
cp $INPUT_FILE $OUTPUT_DIR/USER_INPUT.txt
segfile="$GISTIC_LOC"/run_result/USER_INPUT.txt
refgenefile=$REF_PATH
echo Input file path: $segfile
echo Output directory: $OUTPUT_DIR 
echo Reference file path: $REF_PATH
################


# Run GISTIC
################
echo --- running GISTIC ---
echo Starting gistic...
DOCKER_OUTDIR="$GISTIC_LOC"/run_result
docker run --rm -v $OUTPUT_DIR:$DOCKER_OUTDIR \
  shixiangwang/gistic -b $DOCKER_OUTDIR -seg $segfile -refgene $refgenefile \
  -rx 0 -js 4 -broad 0.98 -cap 1.5 \
  -twoside 1 \
  -genegistic 1 -smallmem 0 -conf 0.99 -armpeel 1 -savegene 1 -saveseg 1 -qvt 0.1
echo Run finished. Hooray!
################