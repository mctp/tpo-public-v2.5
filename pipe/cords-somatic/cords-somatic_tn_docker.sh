#!/usr/bin/env bash
echo $(date) >> /job/docker_start

set -eu
shopt -s nullglob

#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate

source /job/config.txt
PIPECODE="${0%/*}"
NCORES=$(nproc --all)
MEM=$(free -g -h --si |sed -n '2p'|awk '{print $2}'|sed 's/G//')

PFX=/output/$ID

## collect input CRAMs/BAMs
echo $(date) "Collecting input files" >> /job/docker_script.log
BID=$(basename $ALNT)
if [ -f "/input/$BID/$BID.cram" ]; then
    BAMT="/input/$BID/$BID.cram"
elif [ -f "/input/$BID/$BID.bam" ]; then
    BAMT="/input/$BID/$BID.bam"
fi
if [ ! -f "/input/$BID/$BID.cram.crai" ] && [ ! -f "/input/$BID/$BID.bam.bai" ] && \
   [ ! -f "/input/$BID/$BID.bai" ] && [ ! -f "/input/$BID/$BID.crai" ]; then
    $SENTIEON_INSTALL_DIR/bin/sentieon util index $BAMT
fi
BID=$(basename $ALNN)
if [ -f "/input/$BID/$BID.cram" ]; then
    BAMN="/input/$BID/$BID.cram"
elif [ -f "/input/$BID/$BID.bam" ]; then
    BAMN="/input/$BID/$BID.bam"
fi
if [ ! -f "/input/$BID/$BID.cram.crai" ] && [ ! -f "/input/$BID/$BID.bam.bai" ] && \
   [ ! -f "/input/$BID/$BID.bai" ] && [ ! -f "/input/$BID/$BID.crai" ]; then
    $SENTIEON_INSTALL_DIR/bin/sentieon util index $BAMN
fi
BAMT_SAMPLE=$(samtools view -H $BAMT | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
BAMN_SAMPLE=$(samtools view -H $BAMN | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)

## main pipeline
if [ -z "$SENTIEON_LICENSE" ]; then
    IS="-I $BAMT -I $BAMN"
    echo $(date) "Running in opensource mode" >> /job/docker_script.log
    conda activate opensource
    source $PIPECODE/cords-somatic_tn_docker-opensource.sh
else
    IS="-i $BAMT -i $BAMN"
    echo $(date) "Running in sentieon mode" >> /job/docker_script.log
    source $PIPECODE/cords-somatic_tn_docker-sentieon.sh
fi

echo $(date) "Finished." >> /job/docker_script.log
echo $(date) >> /job/docker_done
