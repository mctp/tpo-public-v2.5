#!/bin/bash

# change ID (read group identifier) tag from whatever it was to $PREFIX ($1) and change SM (read group sample) from whatever it was to $SAMPLE ($2) in SAM/BAM/CRAM file 
# old index is removed, new index is generated by default
# usage: fixrg.sh read_group_identifier sample BAM_file

set -eu
PREFIX=$1
SAMPLE=$2
BAM=$3
CORES=$(nproc --all)

BAM_PTH="${BAM%.*}"
TMP_BAM=$BAM_PTH-tmp.bam

# remove old index
if [[ -e ${BAM}.bai ]]; then
    rm ${BAM}.bai
fi
if [[ -e ${BAM}.crai ]]; then
    rm ${BAM}.crai
fi

# replace readgroup in reads, input: cram, output: bam in the pipeline
samtools addreplacerg -r "@RG\tID:$PREFIX\tSM:$SAMPLE\tLB:$SAMPLE" -@$CORES --output-fmt BAM -o $TMP_BAM $BAM
mv $TMP_BAM $BAM

# fix index
samtools index -@$CORES $BAM
