#!/usr/bin/env bash
set -eu

export TPO_ROOT="$(dirname $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ))"

CONFIG=$TPO_ROOT/config/mokefile_grch38.config
REPO=gs://$(grep -Po 'WORK_BUCKET = \K\w(.*)' $CONFIG)/repo

TSAMPLE=$1
NSAMPLE=$2
TALN=$3
NALN=$4
ARGS=${@:5}

## Level 1
$TPO_ROOT/mokefile.py -ll error -config $CONFIG cords_postalign \
                      $TSAMPLE \
                      $TALN &
sleep 1
$TPO_ROOT/mokefile.py -ll error -config $CONFIG cords_postalign \
                      $NSAMPLE \
                      $NALN &
sleep 1
wait
echo "Level 1 done!"

## Level 2
$TPO_ROOT/mokefile.py -ll error -config $CONFIG cords_misc \
                      $TSAMPLE \
                      $REPO/cords-postalign/$TSAMPLE &
sleep 1
$TPO_ROOT/mokefile.py -ll error -config $CONFIG cords_misc \
                      $NSAMPLE \
                      $REPO/cords-postalign/$NSAMPLE &
sleep 1
$TPO_ROOT/mokefile.py -ll error -config $CONFIG cords_somatic \
                      $TSAMPLE.$NSAMPLE \
                      $REPO/cords-postalign/$TSAMPLE \
                      $REPO/cords-postalign/$NSAMPLE &
sleep 1
$TPO_ROOT/mokefile.py -ll error -config $CONFIG cords_structural \
                      $TSAMPLE.$NSAMPLE \
                      $REPO/cords-postalign/$TSAMPLE \
                      $REPO/cords-postalign/$NSAMPLE &
sleep 1
wait
echo "Level 2 done!"

## Level 3
$TPO_ROOT/mokefile.py -ll error -config $CONFIG carat_anno \
    -somatic $REPO/cords-somatic/$TSAMPLE.$NSAMPLE \
    -structural $REPO/cords-structural/$TSAMPLE.$NSAMPLE \
    $TSAMPLE.$NSAMPLE &
sleep 1
$TPO_ROOT/mokefile.py -ll error -config $CONFIG cords_cnvex \
                      -alnt $REPO/cords-postalign/$TSAMPLE \
                      -alnn $REPO/cords-postalign/$NSAMPLE \
                      -tvar $REPO/cords-somatic/$TSAMPLE.$NSAMPLE \
                      $TSAMPLE.$NSAMPLE &
sleep 1
wait
echo "Level 3 done!"
