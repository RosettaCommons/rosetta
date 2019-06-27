#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mp_f19_sequence_recovery/3.plot.py
## @brief this script is part of the franklin2019 sequence recovery test
## @author Sergey Lyskov, Rebecca F. Alford (ralford3@jhu.edu)

import os, sys, subprocess, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#plot metadata
outfile = "plot_results.png"
ncols = 3
nrows = 1
categories = ['all', 'lipid', 'aqueous']

# figure size
width = 7.5 * ncols
height = 6 * nrows

plt.rc("font", size=20)
plt.rcParams['figure.figsize'] = width, height #width, height

# Subplot for Sequence recovery
plt.subplot( nrows, ncols, 1 )
plt.xlabel( "Residue Subset" )
plt.ylabel( "Sequence Recovery (%)" )
plt.title( "Sequence Recovery" )

# scatterplot of the data
plt.bar(np.asarray(categories), np.asarray(recov))
	
# add horizontal and vertical lines for cutoff
plt.axhline(y=float(cutoffs_recov), color='r', linestyle='-')

# Subplot for number non-random
plt.subplot( nrows, ncols, 2 )
plt.xlabel( "Residue Subset" )
plt.ylabel( "Fraction Recovered" )
plt.title( "Non-Random Recovered" )

# scatterplot of the data
plt.bar(np.asarray(categories), np.asarray(non_random_recov))
	
# add horizontal and vertical lines for cutoff
plt.axhline(y=float(cutoffs_non_random), color='r', linestyle='-')

# Subplot for kl-divergence
plt.subplot( nrows, ncols, 3 )
plt.xlabel( "Residue Subset" )
plt.ylabel( "KL Divergence" )
plt.title( "KL-Divergence" )

# scatterplot of the data
plt.bar(np.asarray(categories), np.asarray(kl_divergence))
	
# add horizontal and vertical lines for cutoff
plt.axhline(y=float(cutoffs_kl_div_all), color='r', linestyle='-')
plt.axhline(y=float(cutoffs_kl_div_subset), color='orange', linestyle='-')

#save figure
plt.tight_layout()
plt.savefig( outfile )

benchmark.save_variables('debug targets nstruct working_dir testname results outfile cutoffs_rmsd_dict cutoffs_score_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
