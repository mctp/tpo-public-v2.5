#!/usr/bin/env bash
echo $(date) >> /job/docker_start
echo "Started bcl2fastq." >> /job/docker_script.log

set -eu
#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate

source /job/config.txt
NCORES=$(nproc --all)

# Use custom barcodes if present
if [[ $BCL2FASTQ_SHEET =~ gs://* ]]; then
  BCL2FASTQ_SHEET=/job/input/custom_barcode.csv
fi

if ${REVERSE,,}; then
  Rscript /code/pipe/bcl/reverse_complement.R -i $BCL2FASTQ_SHEET -o /job/input/custom_barcode_rev_comp.csv
  BCL2FASTQ_SHEET=/job/input/custom_barcode_rev_comp.csv
fi

bcl2fastq \
    $BCL2FASTQ_ARGS \
    -l WARNING \
    -p $NCORES \
    -R /flowcell \
    --output-dir=/job/tmp/bcl2fastq \
    --interop-dir=/job/tmp/bcl2fastq/InterOp \
    --stats-dir=/job/tmp/bcl2fastq/Stats \
    --reports-dir=/job/tmp/bcl2fastq/Reports \
    --ignore-missing-positions \
    --ignore-missing-controls \
    --ignore-missing-filter \
    --ignore-missing-bcls \
    --sample-sheet=$BCL2FASTQ_SHEET

rm -rf /input/flowcell

echo "Finished bcl2fastq." >> /job/docker_script.log
