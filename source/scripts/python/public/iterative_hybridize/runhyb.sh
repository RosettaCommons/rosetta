#!/bin/bash
pwd=$1
ROSETTABIN=$2
ROSETTADB=$3
SCRIPTDIR=$4
EXTENSION=$5
iter=$6
templatestr=$7
mode=$8
iseed=$9
wts=${10}
extra=${11}

cd $pwd/iter_$iter

$ROSETTABIN/rosetta_scripts.$EXTENSION -database $ROSETTADB \
    -hb_cen_soft \
    -in:file:fasta ../input.fa \
    -parser:protocol $SCRIPTDIR/$mode.xml -nstruct 1 \
    -out:file:silent_struct_type binary \
    -default_max_cycles 200 \
    -out:file:silent gen.out -out:prefix iter${iter}_from_${mode}.${iseed} \
    -parser:script_vars $templatestr seed=$iseed wts=$wts \
    -fix_disulf ../disulf.def \
    -frag_weight_aligned 0.1 \
    -overwrite \
    $extra \
    >&hybrid.log
