#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# Author: Rahel Frick (frick.rahel@gmail.com)
# Author: Jeliazko Jeliazkov (jeliazkov@jhu.edu)

# TODO Add option for plotting against "ideal" data - JRJ

""" Script for plotting model LHOCs on PDB-derived distribution.

Input:
    - score file output from antibody_H3 app
    - grafted (relaxed) models from antibody.cc app
    - models after H3 refinement (antibody_H3 app)

Output:
    - plot of VL--VH interdomain distance in PDF format
    - plot of VL--VH heavy opening angle in PDF format
    - plot of VL--VH light opening angle in PDF format
    - plot of VL--VH packing angle in PDF format

Example Usage:
  $ python plot_LHOC.py \
        -h3_fasc H3_modeling_scores.fasc \
        -graft_dir ./grafting/ \
        -output_dir ./lhoc_analysis \
        -output_name ab_01

"""
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
import sys
import os
import argparse

from constants import *
from ScoreFile import ScoreFile, LHOCDecoys

names = []
outpath = ''
infiles = []
tempfiles = []
listmode = 0


# argument parsing
def parse_args():
    """ Parse command line options. Default are as follows:
            -h3_fasc H3_modeling_scores.fasc
            -graft_dir ./grafting/
            -output_dir ./lhoc_analysis
            -output_name ab_01
    """

    parser = argparse.ArgumentParser(description=__doc__)

    # TODO ensure h3.fasc is file and exists
    h3_sc_group = parser.add_mutually_exclusive_group()
    h3_sc_group.add_argument('-h3_fasc', default='H3_modeling_scores.fasc', help='What is the h3 score file used for?') # previously -s
    h3_sc_group.add_argument('-h3_fasc_list') # previously -l

    # TODO ensure grafting dir is dir and exists
    graft_group = parser.add_mutually_exclusive_group()
    graft_group.add_argument('-graft_dir', default='grafting', help='Something like script looks for grafted models and calculates their shizz...') # previously -s:temp, should we make it so that the user can specify a basename for the grafted model? Naireeta was having issues with this.
    graft_group.add_argument('-graft_dir_list') # previously -l:temp

    outname_group = parser.add_mutually_exclusive_group()
    outname_group.add_argument('-output_name', default='Antibody 1') # previously -s:names
    outname_group.add_argument('-output_name_list') # previously -l:names

    parser.add_argument('-angles_sc') # previously -angle
    parser.add_argument('-rosetta') # previously -rosetta
    parser.add_argument('-output_dir', default='lhoc_analysis') # previously -out
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # first check if using single files or lists
    if args.h3_fasc_list and args.graft_dir_list and args.output_name_list:
        print 'Detected lists.'
        listmode=1
    # gave less than three lists
    elif args.h3_fasc_list or args.graft_dir_list or args.output_name_list:
        print 'Incompatible input (lists+single files)'
        sys.exit()
    else:
        print 'Detected files.'

    infiles = [args.h3_fasc]
    if args.h3_fasc_list: infiles = args.h3_fasc_list #listmode = 1

    tempfiles = [args.graft_dir]
    if args.graft_dir_list: tempfiles = args.graft_dir_list

    names = [args.output_name]
    if args.output_name_list: names = args.output_name_list

    if args.angles_sc: angles_file = args.angles_sc
    if args.rosetta: rosetta_path = args.rosetta
    if args.output_dir: outpath = args.output_dir


    if listmode:
        with open(infiles, 'r') as f:
            infiles =f.readlines()
        with open(tempfiles, 'r') as f:
            tempfiles =f.readlines()
        with open(names, 'r') as f:
            names =f.readlines()


    for i in range(len(infiles)):
        score_file = ScoreFile(i, infiles, outpath, names).plot_for_all_coordinates(tempfiles, angles_file)
