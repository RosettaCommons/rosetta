#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  ligand_scoring_ranking/3.plot.py
## @brief this script is part of ligand scientific test
## @author Sergey Lyskov
## @author Shannon Smith (shannon.t.smith.1@vanderbilt.edu)

import os, sys, subprocess, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import benchmark
import pandas as pd
import numpy as np

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

x_label = "logKa"
y_label = "score"
infile = 'run_processed_score.csv'
outfile = "plot_results.png"


## get target info from run_processed_score.csv outputted from 2.analyze.py. Plot each target (only 5 points each to see what the outliers look like)
spearman_vals = pd.read_csv('run_Spearman.results', sep='[,,\t, ]+', engine='python')
aa = pd.read_csv(infile, sep='[,,\t, ]+', engine='python')
group = aa.groupby('target')
num = len(aa.drop_duplicates(subset='target')) + 1

#number of subplots--should always be 58 (57 targets + 1 for full dataset)
ncols = 4
nrows = 1
if num < 4:
	ncols = num
else:
	nrows = math.ceil( num / 4 )
# figure size
width = 7.5 * ncols
height = 6 * nrows
plt.rc("font", size=20)
plt.rcParams['figure.figsize'] = width, height #width, height

for i,j in group.__iter__():
	plotdf=group.get_group(i)[['#code','logKa','score']]
	## have to switch scores back to negative values after correlation evaluation
	plotdf['score'] = plotdf['score'].apply(np.negative)

	x = list(plotdf['logKa'])
	y = list(plotdf['score'])
	labels = list(plotdf['#code'])
	for index, line in spearman_vals.iterrows():
		if str(line['#Target']) in labels:
			spearman = line['spearman']

	plt.subplot(nrows, ncols, i )
	plt.xlim(0, round(max(x) + 1))
	plt.ylim(min(y)*1.2, 0)
	plt.xlabel( 'logKa' )
	plt.ylabel( 'score' )
	plt.title( "Target" + ' ' + str(i) )
	plt.annotate("Spearman Correlation = %0.3f"%(spearman), xy=(0,min(y)*1.1), color='blue')
	plt.scatter(x, y)
	for label in range(0,len(x)):
		plt.annotate(str(labels[label]), xy=(x[label], y[label]))

## plot full set. Output Pearson correlation.
aa['score'] = aa['score'].apply(np.negative)
x = list(aa['score'])
y = list(aa['logKa'])
plt.subplot(nrows, ncols, num )
plt.xlim(min(x)*1.2, 0)
plt.ylim(0, round(max(y)+1))
plt.xlabel('score')
plt.ylabel('logKa')
plt.title("Full Dataset")
plt.annotate("Pearson correlation coefficient (R) = %0.3f" % (testr), xy=(min(x)*1.2, 1.2), color='blue')
plt.annotate("Regression: logKa = %.2f + %.2f * score" % ((regr_coef), regr_intercept), xy=(min(x)*1.2,0.2), color='blue')
plt.scatter(x, y)
b = [((regr_coef*i*-1) + (regr_intercept)) for i in x]
plt.plot(x,b, linestyle='--')

#save figure
plt.tight_layout()
plt.savefig( outfile )

benchmark.save_variables('debug targets working_dir testname scorefiles datafile outfile failures below_25 num')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
