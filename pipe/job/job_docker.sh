#!/usr/bin/env bash
echo $(date) >> /job/docker_start

set -eu
shopt -s nullglob

#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate

source /job/config.txt

NCORES=$(nproc --all)

PFX=/output/$ID

chmod +x /job/job_script
/job/job_script /job/job_params

echo $(date) >> /job/docker_done
