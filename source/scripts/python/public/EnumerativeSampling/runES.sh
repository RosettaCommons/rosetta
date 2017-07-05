#!/bin/sh

$ROSETTA3/source/bin/rosetta_scripts.default.linuxgccrelease\
    -database $ROSETTA3/database\
    -parser:protocol $1 \
    -edensity:mapfile $3\
    -s $2 \
    -default_max_cycles 200 \
    -ignore_unrecognized_res \
    -nstruct 1 \
    -overwrite \
    -parser::script_vars readbeams=$4 \
    -parser::script_vars beams=$5 \
    -parser::script_vars steps=$6 \
    -parser::script_vars pcount=$7 \
    -parser::script_vars filterprevious=$8 \
    -parser::script_vars filterbeams=$9 \
    -edensity:sliding_window 3\
    -mapreso 2 \
    -missing_density_to_jump
