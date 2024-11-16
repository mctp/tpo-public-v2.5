#!/usr/bin/env bash
echo $(date) >> /job/docker_start
echo "Started." >> /job/docker_script.log

set -eu
shopt -s nullglob

#### CONDA
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate carat

source /job/config.txt


if [ -z "$MACHINE_NCORES" ]; then
    NCORES=$(nproc --all)
else
    NCORES=$MACHINE_NCORES
fi


for input in /input/*{somatic,germline}-{tinscope,tnscope,mutect2}.vcf.gz; do

    fltvcf=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-flt.vcf/g'))
    incvcf=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-inc.vcf/g'))
    mnvvcf=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-mnv.vcf/g'))
    tmpvcf=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-tmp.vcf/g'))
    vepvcf=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-vep.vcf/g'))
    output=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-annotated.vcf/g'))
    report=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-report/g'))
    
    ## if ANNO_REGIONS is provided restrict variants to those regions
    if [ ! -z "$ANNO_REGIONS" ]; then
        bcftools view $input --threads $NCORES -T $ANNO_REGIONS > $fltvcf 2>> /job/docker_script.log
        bgzip -f -@$NCORES $fltvcf
        tabix $fltvcf.gz
        input=$fltvcf.gz
    fi
    
    ## if ANNO_FILTER_SOMATIC is provided restrict variants which pass this filter
    if [ ! -z "$ANNO_FILTERSOMATIC" ]; then
        bcftools view $input --threads $NCORES $ANNO_FILTERSOMATIC > $incvcf 2>> /job/docker_script.log
        bgzip -f -@$NCORES $incvcf
        tabix $incvcf.gz
        input=$incvcf.gz
    fi
    
    # Phase SNV/Indels into MNVs
    Rscript /code/pipe/carat-anno/carat_mnv.R -g $ANNO_ASSEMBLY -p $NCORES -v $input -o $tmpvcf &>> /job/docker_script.log
    bcftools sort -o $mnvvcf $tmpvcf &>> /job/docker_script.log
    bgzip -f -@$NCORES $mnvvcf
    tabix $mnvvcf.gz

    ## VEP
    vep --dir /tpo/cache/vep --cache --offline --quiet \
        --force_overwrite --vcf --fork $NCORES \
        $ANNO_VEPREFS \
        $ANNO_VEPARGS \
        -i $mnvvcf.gz -o $vepvcf &>> /job/docker_script.log
    bgzip -f -@$NCORES $vepvcf
    tabix -p vcf $vepvcf.gz
    
    ## VCFANNO
    vcfanno -p=$NCORES /code/pipe/carat-anno/$ANNO_CONFIG.toml $vepvcf.gz > $output 2> $output.log.tmp
    grep -v 'chromosome.*not found in' $output.log.tmp > $output.log
    bgzip -f -@$NCORES $output
    tabix $output.gz
     
    Rscript /code/pipe/carat-anno/carat_reportvariants.R \
            -v $output.gz \
            -g $ANNO_ASSEMBLY \
            $ANNO_REPORTARGS \
            -o $report.rds &>> /job/docker_script.log

    ## clean-up
    rm -f $fltvcf*
    rm -f $incvcf*
    rm -f $mnvvcf*
    rm -f $tmpvcf*
    rm -f $vepvcf*
    
    rm -f /output/*vcf_warnings.txt
    rm $output.log.tmp
    
done

for input in /input/*structural-tnscope.vcf.gz; do
    output=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-annotated.vcf/g'))
    tmpvcf=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-tmp.vcf/g'))
    report=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-report/g'))

    ## Add Manta and/or GRIDSS
    mantavcf=$(echo $input | sed 's/tnscope/manta-somaticsv/g')
    gridssvcf=$(echo $input | sed 's/tnscope/gridss-somatic/g')
    if [[ -f "$mantavcf" ]] || [[ -f "$gridssvcf" ]]; then
        Rscript /code/pipe/carat-anno/carat_add_manta_and_gridss.R \
                -t $input \
                -m $mantavcf \
                -d $gridssvcf \
                -g $ANNO_ASSEMBLY \
                -o $tmpvcf &>> /job/docker_script.log
        bgzip -f -@$NCORES $tmpvcf
        tabix $tmpvcf.gz
        input=$tmpvcf.gz
    fi
    
    ## VEP
    vep --dir /tpo/cache/vep --cache --offline --quiet \
        --force_overwrite --vcf --merged --fork $NCORES \
        $ANNO_VEPREFS \
        $ANNO_VEPARGS \
        -i $input -o $output &>> /job/docker_script.log
    bgzip -f -@$NCORES $output
    tabix -p vcf $output.gz
    
    Rscript /code/pipe/carat-anno/carat_reportvariants.R \
            -v $output.gz \
            --structural \
            -g $ANNO_ASSEMBLY \
            $ANNO_REPORTARGS \
            -o $report.rds &>> /job/docker_script.log

    ## clean-up
    rm -f $tmpvcf*
    rm -f /output/*vcf_warnings.txt
    
done


for input in /input/*{somatic,germline}-{dnascope,haplotypecaller}.vcf.gz; do
    
    fltvcf=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-flt.vcf/g'))
    incvcf=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-inc.vcf/g'))
    vepvcf=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-vep.vcf/g'))
    output=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-annotated.vcf/g'))
    report=/output/$(basename $(echo $input | sed 's/\.vcf\.gz$/-report/g'))
    
    ## if ANNO_REGIONS is provided restrict variants to those regions
    if [ ! -z "$ANNO_REGIONS" ]; then
        bcftools view $input --threads $NCORES -T $ANNO_REGIONS > $fltvcf
        bgzip -f -@$NCORES $fltvcf
        tabix $fltvcf.gz
        input=$fltvcf.gz
    fi

    ## if ANNO_FILTERGERMLINE is provided restrict variants which pass this filter
    if [ ! -z "$ANNO_FILTERGERMLINE" ]; then
        bcftools view $input --threads $NCORES $ANNO_FILTERGERMLINE > $incvcf
        bgzip -f -@$NCORES $incvcf
        tabix $incvcf.gz
        input=$incvcf.gz
    fi

    ## VEP
    if [ ! -z "$ANNO_VEPARGS" ]; then
        vep --dir /tpo/cache/vep --cache --offline --quiet \
            --force_overwrite --vcf --merged --fork $NCORES \
            $ANNO_VEPREFS \
            $ANNO_VEPARGS \
            -i $input -o $vepvcf
        bgzip -f -@$NCORES $vepvcf
        tabix -p vcf $vepvcf.gz
        input=$vepvcf.gz
    fi
    
    ## VCFANNO
    vcfanno -p=$NCORES /code/pipe/carat-anno/$ANNO_CONFIG.toml $input > $output 2> $output.log.tmp
    grep -v 'chromosome.*not found in' $output.log.tmp > $output.log
    bgzip -f -@$NCORES $output
    tabix $output.gz
    
    Rscript /code/pipe/carat-anno/carat_reportvariants.R \
            -v $output.gz \
            -g $ANNO_ASSEMBLY \
            $ANNO_REPORTARGS \
            -o $report.rds &>> /job/docker_script.log

    ## clean-up
    rm -f $fltvcf*
    rm -f $incvcf*
    rm -f $vepvcf*
    rm -f /output/*vcf_warnings.txt
    rm $output.log.tmp

done

# done

echo "Finished." >> /job/docker_script.log
echo $(date) >> /job/docker_done
