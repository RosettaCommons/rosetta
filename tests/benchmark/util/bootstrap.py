#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   bootstrap.py
## @brief  Functions for bootstrapping docking results to evaluate quality metrics
## @author Ameya Harmalkar, Shourya S. Roy Burman

import argparse, os, os.path, sys, math, random
import numpy as np
import pandas as pd

#=======================================
def count_high_res_success(sample, n):
    """This functions measures the metric to calculate success. It takes in the list of
    the randomly selected decoys that need to be chosen and returns the success count

    **Arguments**
        sample : *list,list,float,float,str*
        Nested list containing floats and strings .i.e. scores, rms and description
        of the PDB decoys
        
        n : *int*
        The integer of decoys that need to be randomly selected.
        
    **Returns**
        success_count : *int*
        The successful/acceptable decoys present in the random selection from 1000 decoys.
    """
    success_count = 0
    # select CAPRI rank > 0 i.e. acceptable or better quality
    if len(sample) < n :
        if sample[0][1] >= 1:
            success_count += 1
    else:
        for i in range(n):
            if sample[i][1] >= 1:
                success_count += 1
    return success_count


#=======================================
def randomly_sample(scorefile, N_top):
    """Evaluates the <NX> metric where X=N_top, with randomly picked score entries 
    from the scorefile. 

    **Arguments**
        scorefile : *str*
        N_top : *int*
            The top X decoys to be evauated against.
        
    **Returns**
        N_top : *int*
            The <NX> quality metric
    """
    sample = []
    score_file = open(scorefile).readlines()
    headers = score_file[1].split()[1:]
    score_data = [line.split()[1:] for line in score_file[2:]]
    df = pd.DataFrame( score_data, columns=headers, dtype=float )
    
    for i in range(len(score_file) - 2):
        entry = random.randint(1, len(score_file)-2) 
        # appends I_sc, CAPRI_rank
        sample.append( [df.iloc[entry-1]["I_sc"], df.iloc[entry-1]["CAPRI_rank"]])

    sample.sort(key=lambda x:x[0])
    ntop = count_high_res_success(sample, N_top)

    return ntop


#=======================================
def bootstrap_NX(scorefile, N_top, sample_size=1000):
    """To evaluate the average <NX> metric by bootstrapping for defined sample size.
    """

    NX = []

    # sample 1000 times
    for i in range(sample_size):
        ntop = randomly_sample(scorefile, N_top)
        NX.append(ntop)
    
    NX = np.array(NX)
    return np.mean(NX), np.std(NX)
