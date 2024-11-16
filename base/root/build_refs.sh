#!/usr/bin/env bash
set -eu

#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate refs

#### Sentieon
if [ -z "$SENTIEON_LICENSE" ]; then
    BWA=`which bwa`
else
    BWA=$SENTIEON_INSTALL_DIR/bin/bwa
fi

#### CACHE

## update REF_CACHE
mkdir -p /tpo/cache/hts-ref
seq_cache_populate.pl -root /tpo/cache/hts-ref $CORDS_ALIGN_FASTA &
seq_cache_populate.pl -root /tpo/cache/hts-ref $CRISP_ALIGN_FASTA &

## untar VEP_CACHE
mkdir -p /tpo/cache/vep
tar --no-same-owner -xf $CARAT_ANNO_VEPCACHE -C /tpo/cache/vep &

#### DNA

## make BWA indices
mkdir -p /tpo/indices/bwa
cd /tpo/indices/bwa
bedtools maskfasta -fi $CORDS_ALIGN_FASTA -bed $CORDS_ALIGN_MASK -fo $(basename $CORDS_ALIGN_FASTA)
$BWA index -p /tpo/indices/bwa/$CORDS_ALIGN_NAME $(basename $CORDS_ALIGN_FASTA) &
cp $CORDS_ALIGN_FASTA.alt /tpo/indices/bwa/$CORDS_ALIGN_NAME.alt

#### RNA

## make STAR index
mkdir -p /tpo/indices/star/$CRISP_ALIGN_NAME
cd /tpo/indices/star/$CRISP_ALIGN_NAME
bedtools maskfasta -fi $CRISP_ALIGN_FASTA -bed $CRISP_ALIGN_MASK -fo $(basename $CRISP_ALIGN_FASTA)
STAR --runMode genomeGenerate --genomeDir /tpo/indices/star/$CRISP_ALIGN_NAME \
    --genomeFastaFiles $(basename $CRISP_ALIGN_FASTA) \
    --runThreadN $NCORES \
    --sjdbOverhang 125 --sjdbScore 2 --sjdbGTFfile $CRISP_QUANT_GTF &

## make minimap2 index
mkdir -p /tpo/indices/minimap2
minimap2 -G500k -x splice -d /tpo/indices/minimap2/$CRISP_ALIGN_NAME.mmi $CRISP_ALIGN_FASTA &

## make gmap index
mkdir -p /tpo/indices/gmap
gmap_build -d $CRISP_ALIGN_NAME $CRISP_ALIGN_FASTA &

## make kallisto index
mkdir -p /tpo/indices/kallisto
kallisto index -i /tpo/indices/kallisto/$CRISP_QUANT_NAME.idx $CRISP_QUANT_FASTA &

## build msisensor2 index
mkdir -p /tpo/indices/msisensor
msisensor2 scan -d $CORDS_ALIGN_FASTA -o /tpo/indices/msisensor/$CORDS_ALIGN_NAME.msi2 &

## build mi_msi homopolymer GRanges object
mkdir -p /tpo/indices/mi_msi
Rscript /code/base/root/build_msi_ref.R -f $CORDS_ALIGN_FASTA -o /tpo/indices/mi_msi/${CORDS_ALIGN_NAME}_homopolymers.rds &

wait

## fix permissions
chmod 755 -R /tpo/cache/hts-ref