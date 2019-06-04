#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cofactor_binding_sites/3.plot.py
## @brief this script is part of cofactor_binding_sites scientific test
## @author Amanda Loshbaugh
## @author aloshbau@gmail.com

import os, sys, subprocess, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import benchmark
import numpy

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

results_index = {
	'pdb':0,
	'aa':1,
	'design_natural_similarity':2,
	'null_natural_similarity':3,
	'RankTop':4
}

def make_violinplot_data( benchmark ):
	results = open( f"{working_dir}/output/{testname}_results.txt" ).readlines()
	#data_dict = {}

	#data_dict[ benchmark ] = { 'design_data':[], 'null_data':[], 'RankTop':[] }
	data_dict = { 'design_data':[], 'null_data':[], 'RankTop':[] }
	for line in results[1:]:
		
		data_dict['design_data'].append( float(line.split()[ results_index[ 'design_natural_similarity' ] ] ) )
		data_dict['null_data'].append( float(line.split()[ results_index[ 'null_natural_similarity' ] ] ) )
		data_dict['RankTop'].append( float(line.split()[ results_index[ 'RankTop' ] ] ) )

	return data_dict

def format_parts( parts ):
	for pc in parts['bodies']:
		pc.set_edgecolor('black')
	for partname in ('cmins','cmaxes', 'cbars'):
		vp = parts[partname]
		vp.set_edgecolor('black')
		vp.set_linewidth(1)
	for partname in ['cmedians']:
		vp = parts[partname]
		vp.set_edgecolor('black')
		vp.set_linewidth(2)
	return

def violinplot( system, data_for_pdb, fig ):
	failure = False
	
	gs = GridSpec(2,3)
	ax_similarity = fig.add_subplot( gs[0:1,0:2] )
	ax_rank = fig.add_subplot( gs[0:1,2:3] )

	# Profile similarity
	data = [ data_for_pdb['design_data'] ] #, data_for_pdb['null_data'] ]
	parts = ax_similarity.violinplot(data, showmeans=False, showextrema=True, showmedians=True )
	ax_similarity.set_title( 'Profile similarity\n(higher is better)', fontdict={'fontsize':10}, ha='center') #y=1.15, ha='center' )
	ax_similarity.set_xticks( [1,2] )
	ax_similarity.set_xticklabels( [ 'Design' ] ) #[ 'Design', 'Null' ]
	ax_similarity.set_ylabel( 'Similarity to\nNatural Profile' )
	ax_similarity.set_ylim(0,1)
	ax_similarity.set_yticks( numpy.arange(0,1.2,.2) )
	format_parts(parts)
	
	plt.tight_layout()
	
	# Use boxplot to show quartiles
	medianprops = dict(linewidth=0)
	boxprops = dict(linewidth=0)
	bp = ax_similarity.boxplot(data, medianprops=medianprops, boxprops=boxprops, widths=.1, manage_xticks=False, showfliers=False, patch_artist=True )
	for patch in bp['boxes']:
		patch.set_facecolor('black')
	# Median cutoff
	cutoff = 0.7 #numpy.median( data_for_pdb['null_data'] )
	ax_similarity.axhline(y=cutoff, color='g', linestyle=':', linewidth=3, alpha=0.7, label='median cutoff')
	design_median = numpy.median(data_for_pdb['design_data'])
	print( 'design_median = {0}'.format( design_median )  )
	pass_='PASS'
	if design_median < cutoff and config['debug'] == False:
		failure = True
		pass_='FAIL'
	ax_similarity.axhline(y=design_median, color='g', linestyle='-', linewidth=3, alpha=0.7, label='design median ({0})'.format(pass_) )
	# 25th percentile cutoff
	cutoff = 0.5 #numpy.median( data_for_pdb['null_data'] )
	ax_similarity.axhline(y=cutoff, color='b', linestyle=':', linewidth=3, alpha=0.7, label='25th quartile cutoff')
	design_25th_quartile = numpy.percentile( data_for_pdb['design_data'], 25 )
	print( 'design_25th_quartile = {0}'.format( design_25th_quartile ) )
	pass_='PASS'
	if design_25th_quartile < cutoff and config['debug'] == False:
		failure = True
		pass_='FAIL'
	ax_similarity.axhline(y=cutoff, color='b', linestyle='-', linewidth=3, alpha=0.7, label='design 25th quartile ({0})'.format(pass_) )
	# legend
	ax_similarity.legend(loc='right', bbox_to_anchor=(.7, -1), shadow=False, ncol=1, fontsize=12)

	# RankTop
	data = [ data_for_pdb['RankTop'] ]
	parts = ax_rank.violinplot(data, showmeans=False, showextrema=True, showmedians=True )
	ax_rank.set_title( 'Rank Top\n(lower is better)', fontdict={'fontsize':10}, ha='center')
	ax_rank.set_xticks( [1] )
	ax_rank.set_xticklabels( [ 'Design' ] )
	ax_rank.set_ylabel( 'RankTop' )
	ax_rank.set_ylim(0,20)
	ax_rank.set_yticks( numpy.arange(0,21,5) )
	format_parts(parts)
	
	# Use boxplot to show quartiles
	medianprops = dict(linewidth=0)
	boxprops = dict(linewidth=0, color='black')
	bp = ax_rank.boxplot(data, medianprops=medianprops, boxprops=boxprops, widths=.1, manage_xticks=False, showfliers=False, patch_artist=True )
	for patch in bp['boxes']:
		patch.set_facecolor('black')
	# Median cutoff
	cutoff = 2
	ax_rank.axhline(y=cutoff, color='b', linestyle=':', linewidth=3, alpha=0.7, label='median cutoff')
	design_median = numpy.median(data_for_pdb['RankTop'])
	pass_='PASS'
	if design_median > cutoff:
		failure = True
		pass_='FAIL'
	ax_rank.axhline(y=cutoff, color='b', linestyle='-', linewidth=3, alpha=0.7, label='design median ({0})'.format(pass_) )
	# Legend
	ax_rank.legend(loc='right', bbox_to_anchor=(1.1, -1), shadow=False, ncol=1, fontsize=12)
	
	sys_dict = {
		# Cofactor
		'cofactor_3DK9' : 'Glutathione Reductase\nwith cofactor FAD',
		'cofactor_1ZK4' : 'Alcohol Dehydrogenase\nwith cofactor NAP',
		'cofactor_3DLC' : 'Methyltransferase\nwith cofactor SAM',
		'cofactor_2XBN' : 'Aminotransferase\nwith cofactor PMP',
		'cofactor_2IJ2' : 'Cytochrome P450\nwith cofactor HEM',
		'cofactor_3R2Q' : 'Glutathione Transferase\nwith cofactor GSH',
		'cofactor_1F4P' : 'Flavodoxin\nwith cofactor FMN',
		'cofactor_3S6F' : 'Acetyltransferase\nwith cofactor COA',
		'cofactor': 'Cofactor binding sites',
	}
	system = sys_dict[ system ]
	#fig.suptitle('{0} ( n = {1} )'.format( 'Cofactor Binding Sites', len(data[0]) ), y=1.03, ha='center' )#, size=20)

	return failure, len(data[0])


def make_violinplot_pdf():
	benchmark = testname.split('_')[0]

	# initialize figure
	figsize = (6,3) #(7.5,10)
	fig = plt.figure( figsize=figsize )
	#pdf_file = PdfPages( f'{working_dir}/results.pdf' )
	
	# plot
	data_dict = make_violinplot_data( benchmark )
	#for pdb in data_dict:
		#failure = violinplot( pdb, data_dict[pdb], fig )
	failure, n_data_points = violinplot( benchmark, data_dict, fig )
	
	#pdf_file.savefig( bbox_inches='tight' )
	#pdf_file.close()
	
	outfile = "plot_results.png"
	plt.savefig( outfile )
	
	return failure

failure = make_violinplot_pdf()

benchmark.save_variables('targets working_dir testname results failure n_data_points')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
