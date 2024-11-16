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

## main pipeline
if [ -z "$SENTIEON_LICENSE" ]; then
    echo $(date) "Running in opensource mode" >> /job/docker_script.log
    conda activate opensource
    source $PIPECODE/cords-postalign_docker-opensource.sh
else
    echo $(date) "Running in sentieon mode" >> /job/docker_script.log
    source $PIPECODE/cords-postalign_docker-sentieon.sh
fi

echo $(date) "Creating alnids provenance file" >> /job/docker_script.log
for i in "${ALNA[@]}"; do
  ALNID=$(basename $i)
  echo $ALNID >> /output/$ID-alnids
done

echo $(date) "Finished." >> /job/docker_script.log
echo $(date) >> /job/docker_done
