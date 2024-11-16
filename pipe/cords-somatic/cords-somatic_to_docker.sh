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

#determining tumor only or normal only
if [[ "$ALNT" == "None" ]]; then
  BID=$(basename $ALNN)
  MODE='normal'
  echo $(date) "Proceeding with normal sample" >> /job/docker_script.log
else
  BID=$(basename $ALNT)
  MODE='tumor'
  echo $(date) "Proceeding with tumor sample" >> /job/docker_script.log
fi

## collect input CRAMs/BAMs
echo $(date) "Collecting input files" >> /job/docker_script.log
if [ -f "/input/$BID/$BID.cram" ]; then
    BAM="/input/$BID/$BID.cram"
elif [ -f "/input/$BID/$BID.bam" ]; then
    BAM="/input/$BID/$BID.bam"
fi
if [ ! -f "/input/$BID/$BID.cram.crai" ] && [ ! -f "/input/$BID/$BID.bam.bai" ] && \
   [ ! -f "/input/$BID/$BID.bai" ] && [ ! -f "/input/$BID/$BID.crai" ]; then
    $SENTIEON_INSTALL_DIR/bin/sentieon util index $BAM
fi
BAM_SAMPLE=$(samtools view -H $BAM | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)

## main pipeline
if [ -z "$SENTIEON_LICENSE" ]; then
    IS="-I $BAM"
    echo $(date) "Running in opensource mode" >> /job/docker_script.log
    conda activate opensource
    source $PIPECODE/cords-somatic_to_docker-opensource.sh
else
    IS="-i $BAM"
    echo $(date) "Running in sentieon mode" >> /job/docker_script.log
    source $PIPECODE/cords-somatic_to_docker-sentieon.sh
fi

echo $(date) "Finished." >> /job/docker_script.log
echo $(date) >> /job/docker_done
