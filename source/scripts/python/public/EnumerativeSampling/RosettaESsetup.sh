#!/bin/bash
$ROSETTA3/source/bin/grower_prep.default.macosclangrelease \
    -in::file::fasta $2 \
    -pdb $1 \
    -fragsizes 3 9 \
    -fragamounts 100 20 \
