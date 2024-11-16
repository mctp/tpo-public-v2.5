#!/usr/bin/env bash
sudo sh -c '/opt/miniconda3/bin/conda update -n base -c defaults conda && \
            /opt/miniconda3/bin/conda config --add channels bioconda && \
            /opt/miniconda3/bin/conda config --add channels conda-forge && \
            /opt/miniconda3/bin/conda config --set auto_activate_base false'

__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate
