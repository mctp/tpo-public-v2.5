#!/usr/bin/env bash
echo $(date) >> /job/docker_start
echo "Started rename." >> /job/docker_script.log

set -eu
#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate

source /job/config.txt
NCORES=$(nproc --all)

Rscript /code/pipe/bcl/bcl_rename.R -o /job/tmp/bcl2fastq -s /job/input/$(basename $LIB) -f $OUTPUT_FORMAT

mv /job/tmp/bcl2fastq/rename/* /job/output

# record the sizes of the remaining output files
echo "Remaining fastq files:" >> /job/docker_script.log
du -h /tmp/bcl2fastq/*.fastq.gz >> /job/docker_script.log

echo "Finished rename." >> /job/docker_script.log
echo $(date) >> /job/docker_done
