#!/usr/bin/env bash
echo $(date) >> /job/docker_start

set -eu
shopt -s nullglob

#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate crisp

## SETTINGS
source /job/config.txt
NCORES=$(nproc --all)
MEMORY=$(free -m | awk '/^Mem:/{print $2}')
IO_CORES=8
IO_MEMORY=$(( ($MEMORY - 40000) / ($IO_CORES + 2) ))
BB_MEMORY=48000

## IO
OUT=/output
PFX=$OUT/$ID
TMP=$OUT/tmp
mkdir -p $TMP

## Collect Input
echo $(date) "Collecting FASTQ files" >> /job/docker_script.log
IFS=';' read -ra FQ1SA <<< $FQ1S
FQ1SL=""
for i in "${FQ1SA[@]}"; do
    FID1=$(basename $i)
    FQ1SL+=" /input/$FID1"
done
IFS=';' read -ra FQ2SA <<< $FQ2S
FQ2SL=""
for i in "${FQ2SA[@]}"; do
    FID2=$(basename $i)
    FQ2SL+=" /input/$FID2"
done
tmp=$(tempfile -d $TMP) && rm $tmp
fq1cat=$tmp-cat_1.fq.gz
fq2cat=$tmp-cat_2.fq.gz
cat $FQ1SL > $fq1cat &
cat $FQ2SL > $fq2cat &

wait

## BBCUT
tmp=$(tempfile -d $TMP) && rm $tmp
fq1cut=$tmp-cut_1.fq
fq2cut=$tmp-cut_2.fq
bbduk2.sh in=$fq1cat in2=$fq2cat fref=$ALIGN_RRNA \
          stats=/dev/null out=$fq1cut out2=$fq2cut \
          $ALIGN_CUTARGS \
          threads=$NCORES -Xmx"$BB_MEMORY"m overwrite=true &> /dev/null

## KALLISTO
kallisto quant -t $NCORES -i $QUANT_INDEX $QUANT_ARGS -o $TMP $fq1cut $fq2cut
mv $TMP/abundance.tsv $PFX-abundance.tsv
mv $TMP/abundance.h5 $PFX-abundance.h5
mv $TMP/run_info.json $PFX-run_info.json

Rscript /code/pipe/crisp-quant/gene_abundance.R --gtf $QUANT_GTF \
        --inp $PFX-abundance.tsv --out $PFX-gene_abundance.tsv

## clean-up
rm -rf $TMP

echo $(date) >> /job/docker_done
