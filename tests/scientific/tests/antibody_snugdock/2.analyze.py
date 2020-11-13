#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  antibody_snugdock/2.analyze.py
## @brief this script is part of antibody_snugdock scientific test
## @author Jeliazko Jeliazkov

import os, sys, subprocess, math
import numpy as np
import pandas as pd
import benchmark
from benchmark.util import quality_measures as qm

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()
print(config)

results = {}
scorefiles = []
cutoffs_rmsd_dict = {}
cutoffs_discrim_dict = {}
failures = []

# name of output file and cutoff source file
outfile = "result.txt"
cutoffs = "cutoffs"

# scorefiles and logfiles
scorefiles.extend( [ f'{working_dir}/output/{t}/{t}_rmsd.sc' for t in targets ] )

# read cutoffs
cutoff_df = pd.read_csv(cutoffs, sep='\t')
# first col is the number of ** (CAPRI crit.) models in the top10 by I_sc
# second col is the number of * (CAPRI crit.) models in the top10 by I_sc
# final col is the min number of * (CAPRI crit.) models at all

# open results output file
f = open( outfile, "w" )
f.write("target\t" + "top10_#2star\t" + "top10_#1star\t" + "overall_#1star\n") # min rms is of the interface

# go through scorefiles of targets and compare -- very hard coded -- sorry
for i in range( 0, len( scorefiles ) ):
    test_pass = None

    # read in score file, scores are sorted, first one is lowest
    sf = pd.read_csv(scorefiles[i], skiprows=1, sep='\s+')
    
    # sort by I_sc
    sf.sort_values(by=['I_sc'],inplace=True)
    
    # count number of structures in top 10 with capri criteria 2 or better
    nb2 = sum(sf.loc[:,'CAPRI_rank'][:10] >= 2)
    
    # get index for comparison by pdb id
    idx = cutoff_df[cutoff_df.iloc[:,0] == targets[i]].index

    # count number of structures at all with capri criteria 1 or better
    # every structure should pass this
    nb1_all = sum(sf.loc[:,'CAPRI_rank'] >= 1)
    if nb1_all >= int(cutoff_df.loc[idx,'#n1']):
        test_pass = True
    
    # only one of the two criteria below must be passed (limits failures due to noise)

    # compare to expected, if expected is 0, skip -- assumes length 1, but that should be the case
    pass_one = False
    if nb2 > int(cutoff_df.loc[idx,'#top_n2star']) and int(cutoff_df.loc[idx,'#top_n2star']) > 0:
        pass_one = True
    elif int(cutoff_df.loc[idx,'#top_n2star']) <= 0:
        pass_one = True
    
    # count number of structures in top 10 with capri criteria 1 or better
    nb1 = sum(sf.loc[:,'CAPRI_rank'][:10] >= 1)
    # compare to expected, if expected is 0, skip -- assumes length 1 again
    pass_two = False
    if nb1 > int(cutoff_df.loc[idx,'#top_n1star']) and int(cutoff_df.loc[idx,'#top_n1star']) > 0:
        pass_two = True
    elif int(cutoff_df.loc[idx,'#top_n1star']) <= 0:
        pass_two = True
   
    if not (pass_one or pass_two):
        test_pass= False

    if not test_pass:
        failures.append( targets[i] )
    
    f.write(f'{targets[i]}\t{nb2}\t{nb1}\t{nb1_all}\n')
    target_results = {'nb2':nb2, 'nb1':nb1, 'nb1_all':nb1_all}
    results.update( {targets[i] : target_results} )

f.close()

benchmark.save_variables('debug targets nstruct working_dir testname results scorefiles cutoffs_discrim_dict cutoffs_rmsd_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
