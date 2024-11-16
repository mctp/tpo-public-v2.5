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

## check and correct quality encoding
FQ1L=/input/$(basename $FQ1)
FQ2L=/input/$(basename $FQ2)
if [ -z "$ALIGN_QUAL" ]; then
    ALIGN_QUAL=$(testformat.sh $FQ1L $FQ2L | cut -f1 | uniq)
fi
if [ "$ALIGN_QUAL" != "sanger" ]; then
    echo $(date) "Running Reformatting" >> /job/docker_script.log
    reformat.sh in=$FQ1L in2=$FQ2L out=/tmp/1.fq out2=/tmp/2.fq qin=auto qout=sanger t=$NCORES &>> /job/docker_script.log
    FQ1L=/tmp/1.fq
    FQ2L=/tmp/2.fq
fi

## trim adapters from reads
if [ "$ALIGN_TRIM" != "none" ]; then
    echo $(date) "Running Trimming" >> /job/docker_script.log
    bbduk.sh -Xmx1g t=$NCORES \
        in1=$FQ1L in2=$FQ2L out1=/tmp/1t.fq out2=/tmp/2t.fq ref=/opt/bbmap/resources/$ALIGN_TRIM.fa.gz \
        ktrim=r k=23 mink=11 hdist=1 ignorebadquality=t qin=33 \
        tpe tbo &>> /job/docker_script.log
    rm -f /tmp/1.fq /tmp/2.fq
    FQ1L=/tmp/1t.fq
    FQ2L=/tmp/2t.fq
fi

## main pipeline
if [ -z "$SENTIEON_LICENSE" ]; then
    echo $(date) "Running in opensource mode" >> /job/docker_script.log
    conda activate opensource
    source $PIPECODE/cords-align_docker-opensource.sh
else
    echo $(date) "Running in sentieon mode" >> /job/docker_script.log
    source $PIPECODE/cords-align_docker-sentieon.sh
fi

echo $(date) "Finished." >> /job/docker_script.log
echo $(date) >> /job/docker_done
