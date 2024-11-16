#!/usr/bin/env bash
echo $(date) >> /job/docker_start

set -eu
shopt -s nullglob

#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate crisp
source /code/pipe/crisp-align/crisp-align_helper.sh

source /job/config.txt

NCORES=$(nproc --all)

mkdir -p /input/$ID

PFX=/output/$ID

merge-crisp-align $ALN /input /input/$ID/$ID $ALIGN_FASTA
Rscript /code/pipe/crisp-codac/codac_anno.R -g hg38 -j 4 $CODAC_GTF /input/annotation.rds
Rscript /code/pipe/crisp-codac/codac_run.R -a /input/annotation.rds -c $CODAC_CONFIG /input/$ID $PFX
Rscript /code/pipe/crisp-codac/codac_stat.R \
        -a /input/annotation.rds -c $CODAC_CONFIG $PFX-run.rds $PFX
Rscript /code/pipe/crisp-codac/codac_callsv.R \
        $CODAC_INDEX \
        -a /input/annotation.rds -c $CODAC_CONFIG $PFX-run.rds $PFX

echo $(date) >> /job/docker_done
