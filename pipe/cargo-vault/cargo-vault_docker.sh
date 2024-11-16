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

if [ -z "$CASE" ]; then
    CASE_INPUT=""
else
    CASE_INPUT="--case $CASE"
fi
if [ -z "$COHORT" ]; then
    COHORT_INPUT=""
else
    COHORT_INPUT="--cohort $COHORT"
fi
if [ -z "$ID" ]; then
    ID_INPUT=""
else
    ID_INPUT="--id $ID"
fi
if [ -z "$CARAT_ANNO" ]; then
    CARAT_INPUT=""
else
    CARAT_INPUT="--carat /input/carat-anno"
fi
if [ -z "$CORDS_MISC_TUMOR" ]; then
    MISC_INPUT_TUMOR=""
else
    MISC_INPUT_TUMOR="--misc_t /input/cords-misc-tumor"
fi
if [ -z "$CORDS_MISC_NORMAL" ]; then
    MISC_INPUT_NORMAL=""
else
    MISC_INPUT_NORMAL="--misc_n /input/cords-misc-normal"
fi
if [ -z "$CORDS_CNVEX" ]; then
    CNVEX_INPUT=""
else
    CNVEX_INPUT="--cnvex /input/cords-cnvex"
fi
if [ -z "$CRISP_QUASR_TUMOR" ]; then
    QUASR_INPUT_TUMOR=""
else
    QUASR_INPUT_TUMOR="--tquasr /input/crisp-tquasr"
fi
if [ -z "$CRISP_QUASR_NORMAL" ]; then
    QUASR_INPUT_NORMAL=""
else
    QUASR_INPUT_NORMAL="--nquasr /input/crisp-nquasr"
fi
if [ -z "$CRISP_CODAC" ]; then
    CODAC_INPUT=""
else
    CODAC_INPUT="--codac /input/crisp-codac"
fi
if [ "$VAULT_DUMP" = true ]; then
    DUMP_OUTPUT="--dump $PFX-vault-dump.rds"
else
    DUMP_OUTPUT=""
fi
if [ -z "$VAULT_HOMOPOLYMERS" ]; then
    HOMO_INPUT=""
else
    HOMO_INPUT="--homo $VAULT_HOMOPOLYMERS"
fi
if [ -z "$VAULT_TARGETS" ]; then
    TARGETS_INPUT=""
else
    TARGETS_INPUT="--targets $VAULT_TARGETS"
fi

Rscript /code/pipe/cargo-vault/cargo_vault.R \
        -o $PFX-vault.rds \
        -g $VAULT_GENEGTF \
        -s $VAULT_SETTINGS \
        $HOMO_INPUT \
        $TARGETS_INPUT \
        $DUMP_OUTPUT \
        $CASE_INPUT \
        $COHORT_INPUT \
        $ID_INPUT \
        $CARAT_INPUT \
        $MISC_INPUT_TUMOR \
        $MISC_INPUT_NORMAL \
        $CNVEX_INPUT \
        $QUASR_INPUT_TUMOR \
        $QUASR_INPUT_NORMAL \
        $CODAC_INPUT

echo $(date) >> /job/docker_done
