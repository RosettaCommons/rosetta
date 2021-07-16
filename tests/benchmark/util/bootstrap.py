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
## @edits Morgan Nance (@mlnance; for generality of bootstrap calculations)

import argparse, os, os.path, sys, math, random
import numpy as np
import pandas as pd
import benchmark
from benchmark.util import scorefile_io_utils as sciu

#=======================================
def count_success(sub_sample_df, score, N_top, metric, cutoff, gt_eq_cutoff):
    """This functions measures the metric to calculate success. It takes in the list of
    the randomly selected decoys that need to be chosen and returns the success count

    **Arguments**
        sub_sample_df : *DataFrame*
            A DataFrame containing 'score' and 'metric' for decoy of the random sub-selection.
        
        score : *str*
            The name of the score on which to sort the decoys.

        N_top : *int*
            The top X decoys to be evauated against.

        metric : *str*
            The name of the metric with which to determine docking success.

        cutoff : *float*
            The float cutoff value defining docking success
        
        gt_eq_cutoff : *bool*
            Is success defined by the metric being >= to the cutoff (True)? Or
            is it defined by the metric being <= to the cutoff (False)?
        
    **Returns**
        success_count : *int*
            The successful decoys present in the random sub-selection of decoys.
    """
    # example success cutoff of 1 when gt_eq_cutoff is True:
    # --select CAPRI_rank >= 1 i.e. acceptable or better quality
    # example success cutoff of 2.0 when gt_eq_cutoff is False:
    # --select rmsd <= 2.0 i.e. a near-native decoy
    sub_sample_df = sub_sample_df.sort_values(score)
    if gt_eq_cutoff:
        return sum(sub_sample_df.iloc[:N_top][metric] >= cutoff)
    else:
        return sum(sub_sample_df.iloc[:N_top][metric] <= cutoff)


#=======================================
def randomly_sample(df, sample_size, score, N_top, metric, cutoff, gt_eq_cutoff):
    """Evaluates the <NX> metric where X=N_top, with randomly picked score entries 
    from the scorefile. 

    **Arguments**
        df : *DataFrame*
            The scorefile turned into a Pandas DataFrame

        sample_size : *int*
            The number of random models to select (with replacement) from the original set.
            Ex: len(df) == 10000 & sample_size == 1000, select 1000 random models from the 10000.
                where "with replacement" means you could select the same model more than once.

        score : *str*
            The name of the score on which to sort the decoys.

        N_top : *int*
            The top X decoys to be evauated against.

        metric : *str*
            The name of the metric with which to determine docking success.

        cutoff : *float*
            The float cutoff value defining docking success
        
        gt_eq_cutoff : *bool*
            Is success defined by the metric being >= to the cutoff (True)? Or
            is it defined by the metric being <= to the cutoff (False)?
        
    **Returns**
        N_top : *int*
            The <NX> quality metric
    """
    # examples:
    # --samples "I_sc" (score) and its "CAPRI_rank" (metric)
    # --samples "interaction_energy" (score) and its "rmsd" (metric)

    if len(df) < sample_size:
        # should only be here in debug mode, really
        # debug mode nstruct is very low, say 5-20
        # so should not take 1000 samples if len(df) is low
        sample_size = len(df)

    # Note: could set a random_state for pandas' random sample
    sub_sample_df = df.sample(n=sample_size, replace=True)

    return count_success(sub_sample_df, score, N_top, metric, cutoff, gt_eq_cutoff)


#=======================================
def bootstrap_NX(scorefile,
                 score_str, N_top,
                 metric_str, cutoff, gt_eq_cutoff,
                 sample_size=1000, # N random models selected form orig model set
                 n_samples=5000): # B in doi.org/10.1371
    """To evaluate the average <NX> metric by bootstrapping for defined sample size.
    """
    df = sciu.scorefile_to_dataframe(scorefile)
    # only care about the score_str and metric_str cols in df from scorefile
    df = df[[score_str, metric_str]]

    # evaluate bootstrap <NX> by taking    
    # <n_samples> random, independent sub-samples of size <sample_size>
    NX = []
    for ii in range(n_samples):
        ntop = randomly_sample(df, sample_size,
                               score_str, N_top, # top-N ranked by score
                               metric_str, cutoff, gt_eq_cutoff) # is metric >= or <= cutoff
        NX.append(ntop)
    
    NX = np.array(NX)
    return np.mean(NX), np.std(NX)
