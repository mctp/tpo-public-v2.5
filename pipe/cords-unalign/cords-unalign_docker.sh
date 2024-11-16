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

FILE=$(basename $ALN)
EXT="${FILE##*.}"
BAM=/input/$FILE
echo $(date) "Running bamtofastq" >> /job/docker_script.log
bamtofastq \
    collate=1 \
    exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
    filename=$BAM \
    gz=0 \
    inputformat=$EXT \
    reference=$ALIGN_FASTA \
    T=/output/$ID-tmp \
    F="$PFX"_1.fq \
    F2="$PFX"_2.fq \
    S="$PFX"_s.fq \
    O="$PFX"_o1.fq \
    O2="$PFX"_o2.fq \
    tryoq=1 &>> /job/docker_script.log

ls -lh /output/*.fq >> /job/docker_script.log
rm "$PFX"_s.fq "$PFX"_o1.fq "$PFX"_o2.fq
pigz -p $NCORES "$PFX"_1.fq
pigz -p $NCORES "$PFX"_2.fq

echo $(date) "Finished." >> /job/docker_script.log
echo $(date) >> /job/docker_done
