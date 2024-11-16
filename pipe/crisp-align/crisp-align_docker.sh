#!/usr/bin/env bash
echo $(date) >> /job/docker_start

set -eu
shopt -s nullglob

#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate crisp

source /job/config.txt
NCORES=$(nproc --all)
MEMORY=$(free -m | awk '/^Mem:/{print $2}')
IO_CORES=8
BB_MEMORY=48000

OUT=/output
PFX=$OUT/$ID
TMP=$OUT/tmp
mkdir -p $TMP

#### create input files for STAR
FQ1L=/input/$(basename $FQ1)
FQ2L=/input/$(basename $FQ2)
tmp=$(tempfile -d $TMP) && rm $tmp

## UMI
if [ ! -z "$ALIGN_UMI" ]; then
    fq1umi=$tmp-umi_1.fq
    fq2umi=$tmp-umi_2.fq
    $SENTIEON_INSTALL_DIR/bin/sentieon umi extract \
        $ALIGN_UMI \
        $FQ1L $FQ2L | \
        paste -d $'\035' - - - - - - - - | tee \
        >(cut -d $'\035' -f 1-4 | tr $'\035' "\n" | cut -f 1-2 | tr "\t" "_" > $fq1umi) | \
          cut -d $'\035' -f 5-8 | tr $'\035' "\n" | cut -f 1-2 | tr "\t" "_" > $fq2umi
    FQ1L=$fq1umi
    FQ2L=$fq2umi
fi

## BBCUT
fq1cut=$tmp-cut_1.fq
fq2cut=$tmp-cut_2.fq
cutlog=$PFX-cut.log
cutout=$PFX-cut.out
bbduk2.sh in=$FQ1L in2=$FQ2L fref=$ALIGN_RRNA \
          stats=$cutlog out=$fq1cut out2=$fq2cut \
          $ALIGN_CUTARGS \
          threads=$NCORES -Xmx"$BB_MEMORY"m overwrite=true &> $cutout

## BBMERGE
fq1mrg=$tmp-mrg_1.fq
fq2mrg=$tmp-mrg_2.fq
fq3mrg=$tmp-mrg_3.fq
mrglog=$PFX-mrg.log
mrgout=$PFX-mrg.out
bbmerge.sh in1=$fq2cut in2=$fq1cut out=$fq3mrg \
           outu1=$fq2mrg outu2=$fq1mrg ihist=$mrglog \
           $ALIGN_MRGARGS usejni=t \
           threads=$NCORES -Xmx"$BB_MEMORY"m overwrite=true &> $mrgout

#### STAR alignment

## linear STAR alignment
mkdir -p $TMP/star-alig $TMP/sort
cd $TMP/star-alig
( STAR \
    --readFilesIn $fq2cut $fq1cut \
    --genomeDir $ALIGN_INDEX \
    --genomeLoad "LoadAndKeep" \
    --outSAMattrRGline ID:$ID SM:$SAMPLE LB:$SAMPLE PL:ILLUMINA \
    --runThreadN $NCORES \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 3 \
    --scoreGenomicLengthLog2scale 0 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outStd "SAM" \
    --outReadsUnmapped "Fastx" \
    --outFilterType "BySJout" 2>> /job/docker_script.log ) \
    | $SENTIEON_INSTALL_DIR/bin/sentieon util sort \
                                         -r $ALIGN_FASTA -t $NCORES --temp_dir $TMP/sort \
                                         --cram_write_options version=3.0,compressor=rans \
                                         --sam2bam -o $PFX-alig.cram -i - &>> /job/docker_script.log
pigz -p $((NCORES / 2)) Unmapped.out.mate1 &
pigz -p $((NCORES / 2)) Unmapped.out.mate2 &
wait
mv Unmapped.out.mate1.gz $PFX-unmapped_1.fq.gz
mv Unmapped.out.mate2.gz $PFX-unmapped_2.fq.gz
gzip -c SJ.out.tab > $PFX-sj.tab.gz
mv Log.final.out $PFX-alig.log
cd $OUT

## chimeric STAR alignment PE
mkdir -p $TMP/star-chim-pe
cd $TMP/star-chim-pe
STAR \
     --readFilesIn $fq2mrg $fq1mrg \
     --genomeDir $ALIGN_INDEX \
     --genomeLoad "LoadAndKeep" \
     --runThreadN $NCORES \
     --outStd "SAM" \
     --outReadsUnmapped "None" \
     --alignTranscriptsPerReadNmax 100000 \
     --outSAMattrRGline ID:$ID SM:$SAMPLE LB:$SAMPLE PL:ILLUMINA \
     --outFilterType Normal \
     --alignIntronMax 150000 \
     --alignMatesGapMax 150000 \
     --chimSegmentMin 10 \
     --chimJunctionOverhangMin 1 \
     --chimScoreSeparation 0 \
     --chimScoreJunctionNonGTAG 0 \
     --chimScoreDropMax 1000 \
     --chimScoreMin 1 \
     > /dev/null 2> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon util sort \
                                   -r $ALIGN_FASTA -t $NCORES --temp_dir $TMP/sort \
                                   --cram_write_options version=3.0,compressor=rans \
                                   --sam2bam -o $PFX-chim-pe.cram -i Chimeric.out.sam &>> /job/docker_script.log
mv Log.final.out $PFX-chim-pe.log
gzip -c Chimeric.out.junction > $PFX-chim-pe.jnc.gz
cd $OUT

## chimeric STAR alignment SE
mkdir -p $TMP/star-chim-se
cd $TMP/star-chim-se
STAR \
     --readFilesIn $fq3mrg \
     --genomeDir $ALIGN_INDEX \
     --genomeLoad "LoadAndRemove" \
     --runThreadN $NCORES \
     --outStd SAM \
     --outReadsUnmapped None \
     --alignTranscriptsPerReadNmax 100000 \
     --outSAMattrRGline ID:$ID SM:$SAMPLE LB:$SAMPLE PL:ILLUMINA \
     --outFilterType Normal \
     --alignIntronMax 150000 \
     --alignMatesGapMax 150000 \
     --chimSegmentMin 10 \
     --chimJunctionOverhangMin 1 \
     --chimScoreSeparation 0 \
     --chimScoreJunctionNonGTAG 0 \
     --chimScoreDropMax 1000 \
     --chimScoreMin 1 \
     > /dev/null 2> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon util sort \
                                         -r $ALIGN_FASTA -t $NCORES --temp_dir $TMP/sort \
                                         --cram_write_options version=3.0,compressor=rans \
                                         --sam2bam -o $PFX-chim-se.cram -i Chimeric.out.sam &>> /job/docker_script.log
mv Log.final.out $PFX-chim-se.log
gzip -c Chimeric.out.junction > $PFX-chim-se.jnc.gz
cd $OUT

## Genotyping
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $PFX-alig.cram \
                                   --algo RNASplitReadsAtJunction --reassign_mapq 255:60 $TMP/split.bam &>> /job/docker_script.log
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $ALIGN_FASTA -t $NCORES -i $TMP/split.bam \
                                   --interval $ALIGN_GENOTYPE \
                                   --algo DNAscope --trim_soft_clip --call_conf 20 --emit_conf 20 --emit_mode gvcf \
                                   $PFX-genotype.gvcf &>> /job/docker_script.log
Rscript /code/pipe/common/genotype_calling.R -i $PFX-genotype.gvcf -o $PFX-genotype.csv &>> /job/docker_script.log

## clean-up
rm -rf $TMP

####
echo $(date) >> /job/docker_done
