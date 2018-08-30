#!/bin/sh

#$ROSETTA3/source/bin/rosetta_scripts.default.linuxgccrelease \
ROSETTA/dgdp/source/cmake/build_release/rosetta_scripts \
    -database $ROSETTA3/database \
    -parser:protocol RosettaES.xml \
    -s 1rpb_symm_INPUT.pdb \
    -default_max_cycles 200 \
    -ignore_unrecognized_res \
    -nstruct 1 \
    -overwrite \
    -missing_density_to_jump
