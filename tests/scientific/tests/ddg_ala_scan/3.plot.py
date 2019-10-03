#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/3.plot.py
## @brief this script is part of cartesian_relax scientific test
## @author Sergey Lyskov

import os, sys, subprocess, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import benchmark

from collections import namedtuple
LinregressResult = namedtuple('LinregressResult', ('slope', 'intercept',
                                                   'rvalue', 'pvalue',
                                                   'stderr'))

results={}
#debug targets nstruct working_dir testname results  failures
benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals


config = benchmark.config()

# inputs are header labels from the scorefile to plot, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
outfile = f"{working_dir}/plot_results.png"

score_function = 'talaris2014'
ddg_df = pd.read_csv(f"{working_dir}/output/calc_vs_exp.csv")
ddg_df.plot('DDG_obs', 'DDG_calc', kind='scatter')
fit = LinregressResult(*results[score_function]['fit'])

R   = results[score_function]['R']
MAE = results[score_function]['MAE']
FCC = results[score_function]['FCC']

# Plot data points
y_pred = fit.intercept + fit.slope*ddg_df.DDG_obs
plt.plot(ddg_df.DDG_obs,y_pred, color="black", label="Fitted line")

plt.xlabel('Experimental ddG')
plt.ylabel('Calculated ddG')
plt.suptitle(f'{score_function}: R = {R:2.2f}; MAE = {MAE:2.2f}; FCC = {FCC:2.2f}')

#save figure
#plt.tight_layout()
plt.savefig( outfile )

benchmark.save_variables('debug targets nstruct working_dir testname results outfile cutoffs_rmsd_dict cutoffs_score_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
