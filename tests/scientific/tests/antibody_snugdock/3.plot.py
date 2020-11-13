#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  antibody_snugdock/3.plot.py
## @brief SnugDock tests for CAPRI quality of models in the top 10 (by Isc), plot this data
## @author Sergey Lyskov, modified by Jeliazko Jeliazkov

import os, sys, subprocess, math
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

cc_qual_fig = "cc_qual.png"
funnel_fig = "funnel.png"

# previous output file is results.txt
results = pd.read_csv("result.txt", sep='\t')
cutoffs = pd.read_csv("cutoffs", sep='\t')

# creates an ugly dataframe with bad labels, but ...
# cols starting with '#' are the cutoffs, other cols are the results
merge = cutoffs.merge(results, left_on='#pdb', right_on='target')
merge.index = merge.target

# want plots of three comparisons

# figure size
ncols = 3
nrows = 1
width = 7.5 * ncols
height = 6 * nrows

plt.rc("font", size=20)
plt.rcParams['figure.figsize'] = width, height #width, height
fig = plt.figure()

ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

# number of models with capri quality >= 2 in top 10 vs. cutoff
sdf2 = merge.loc[:,['target','#top_n2star', 'top10_#2star']]
sdf2.rename(columns={'#top_n2star':'cutoff', 'top10_#2star':'result'},inplace=True)
sdf2.plot.bar(ax=ax1)
ax1.set_ylabel("count")
ax1.set_title(">Medium CAPRI Quality Models in (Top 10)")

# number of models with capri quality >= 1 in top 10 vs. cutoff
sdf1 = merge.loc[:,['target','#top_n1star', 'top10_#1star']]
sdf1.rename(columns={'#top_n1star':'cutoff', 'top10_#1star':'result'},inplace=True)
sdf1.plot.bar(ax=ax2)
ax2.set_ylabel("count")
ax2.set_title(">Acceptable CAPRI Quality Models (Top 10)")

# number of models with capri quality >= 1 in all vs. cutoff
sdfall = merge.loc[:,['target','#n1', 'overall_#1star']]
sdfall.rename(columns={'#n1':'cutoff', 'overall_#1star':'result'},inplace=True)
sdfall.plot.bar(ax=ax3)
ax3.set_ylabel("count")
ax3.set_title(">Acceptable CAPRI Quality Models (All)")

plt.tight_layout()
fig.savefig(cc_qual_fig)

# now for the funnel plots
#number of subplots, dictated by columns (rows are flexible)
ncols = 3
nrows = 1
if len( targets ) < ncols:
    ncols = len( targets )
else:
    nrows = math.ceil( len( targets ) / ncols )

# figure size
width = 7.5 * ncols
height = 6 * nrows

plt.rc("font", size=20)
plt.rcParams['figure.figsize'] = width, height #width, height
plt.figure()

# go through scorefiles
for i in range( 0, len( scorefiles ) ):

    # read in score file
    sf = pd.read_csv(scorefiles[i], skiprows=1, sep='\s+')
    
    # create subplot
    plt.subplot( nrows, ncols, i+1 )
    
    # x and y labels
    plt.xlabel( 'Interface RMSD (A)' )
    plt.ylabel( "Interface Score (REU)" )
    
    # set title
    plt.title( targets[i] )
    
    # scatterplot of the data
    plt.plot(sf['Irms'], sf['I_sc'], 'ko')

#save figure
plt.tight_layout()
plt.savefig( funnel_fig )

benchmark.save_variables('debug targets nstruct working_dir testname results outfile cutoffs_rmsd_dict cutoffs_discrim_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
