#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  ligand_scoring_ranking/2.analyze.py
## @brief this script is part of ligand scientific test
## @author Sergey Lyskov
## @author Shannon Smith (shannon.t.smith.1@vanderbilt.edu)

import os, sys, subprocess, math
import numpy as np
import pandas as pd
import benchmark
from benchmark.util import quality_measures as qm
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
import scipy

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()

scorefiles = []
#logfiles = []
failures = []
successes = []

# inputs are header labels from the scorefile, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
x_label = "ligand_rms_no_super_X"
y_label = "interface_delta_X"
datafile = "datafile.txt"
outfile = "result.txt"

# scorefiles and logfiles
scorefiles.extend( [ f'{working_dir}/output/{t}/{t}.score' for t in targets ] )
#logfiles.extend( [ f'{working_dir}/hpc-logs/hpc.{testname}-{t}.*.log' for t in targets ] )

# get column numbers from labels, 1-indexed
x_index = str( subprocess.getoutput( "grep " + x_label + " " + scorefiles[0] ).split().index( x_label ) + 1 )
y_index = str( subprocess.getoutput( "grep " + y_label + " " + scorefiles[0] ).split().index( y_label ) + 1 )

# open results output file
with open( datafile, "w" ) as f:
	f.write("#code score\n")

	# go through scorefiles of targets
	for i in range( 0, len( scorefiles ) ):

		target_results = {}

		# read in score file, scores are sorted, first one is lowest
		x = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | sort -nk2 | awk '{print $" + x_index + "}'" ).splitlines()
		y = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | sort -nk2 | awk '{print $" + y_index + "}'" ).splitlines()
		# map values to floats (were strings)
		#x = list( map( float, x ))
		y = sorted(list( map( float, y )))[0]
		f.write(targets[i] + ' ' + str(y) + '\n')
f.close()
# Now run correlation tests with Kd
with open('CoreSet.dat','r') as f:
	with open('cstemp','w+') as f1:
		for i in f.readlines():
			if i.startswith('#'):
				if i.startswith('#code'):
					f1.writelines(i)
				else:
					continue
			else:
				f1.writelines(i)
	f1.close()
f.close()
aa = pd.read_csv('cstemp', sep='[,,\t, ]+', engine='python')
aa = aa.drop_duplicates(subset=['#code'],keep='first')
bb = pd.read_csv(datafile, sep='[,,\t, ]+', engine='python')
testdf1 = pd.merge(aa,bb,on='#code')
testdf1['score'] = testdf1['score'].apply(np.negative)
testdf2 = testdf1[testdf1.score > 0]
testdf2.to_csv('run_processed_score.csv', columns=['#code', 'logKa', 'score', 'target'], sep='\t', index=False)

regr = linear_model.LinearRegression()
regr.fit(testdf2[['score']], testdf2[['logKa']])
testpredy = regr.predict(testdf2[['score']])
testr = scipy.stats.pearsonr(testdf2['logKa'].values, testdf2['score'].values)[0]
testmse = mean_squared_error(testdf2[['logKa']], testpredy)
num= testdf2.shape[0]
testsd = np.sqrt((testmse*num)/(num-1))
regr_coef = float(regr.coef_)
regr_intercept = float(regr.intercept_)
testr = float(testr)

if os.path.exists('cstemp'):
	os.remove('cstemp')

## done with scoring test
## Ranking test

#Get the representative complex in each cluster
group=testdf1.groupby('target')

def top(df,n=1,column='logKa'):
	return df.sort_values(by=column)[-n:]
toptardf=testdf1.groupby('target').apply(top)
targetlst=toptardf['#code'].tolist()

#Calculate the Spearman correlation coefficient
spearman=pd.DataFrame(index=targetlst,columns=['spearman'])
rankresults=pd.DataFrame(index=range(1,len(targetlst)+1),columns=['Target','Rank1','Rank2','Rank3','Rank4','Rank5'])
tmp=1
for i,j in group.__iter__():
	testdf2=group.get_group(i)[['#code','logKa','score', 'target']]
	testdf2=testdf2.sort_values('score',ascending=False)
	tartemp=top(testdf2)['#code'].tolist()
	tar=''.join(tartemp)
	if len(testdf2) == 5:
		spearman.ix[tar]['spearman']=testdf2.corr('spearman')['logKa']['score']
		rankresults.ix[tmp]['Rank1']=''.join(testdf2[0:1]['#code'].tolist())
		rankresults.ix[tmp]['Rank2']=''.join(testdf2[1:2]['#code'].tolist())
		rankresults.ix[tmp]['Rank3']=''.join(testdf2[2:3]['#code'].tolist())
		rankresults.ix[tmp]['Rank4']=''.join(testdf2[3:4]['#code'].tolist())
		rankresults.ix[tmp]['Rank5']=''.join(testdf2[4:5]['#code'].tolist())
		rankresults.ix[tmp]['Target']=tar
		tmp+=1
	else:
		spearman.drop(tar,inplace=True)

#Print the output of ranking power evluation
spearmanmean=float(spearman['spearman'].sum()/float(spearman.shape[0]))
spearman.to_csv('run_Spearman.results',sep='\t',index_label='#Target')

#Output results
with open(outfile, 'w+') as f2:
	f2.writelines("Summary of the scoring power: ===================================")
	f2.writelines("\nThe regression equation: logKa = %.2f + %.2f * Score" % ((regr.coef_), regr.intercept_))
	f2.writelines("\nNumber of favorable sample (N) = %d"%(num))
	f2.writelines("\nPearson correlation coefficient (R) = %0.3f"%(testr))
	f2.writelines("\nStandard deviation in fitting (SD) = %0.2f"%(testsd))
	f2.writelines("\n=================================================================")
	f2.writelines("\nSummary of the ranking power: ===========================================")
	f2.writelines("\nThe Spearman correlation coefficient (SP) = %0.3f"%(spearmanmean))
f2.close()

if os.path.exists('cstemp'):
	os.remove('cstemp')

## Find failure clusters for ranking test
cc = pd.read_csv('run_Spearman.results', delim_whitespace=True, engine='python')
for index, row in cc.iterrows():
	if float(row['spearman']) < 0.25:
		failures.append(row['#Target'])
	else:
		successes.append(row['#Target'])

benchmark.save_variables('debug targets working_dir testname scorefiles datafile outfile failures successes regr_coef regr_intercept testr')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
