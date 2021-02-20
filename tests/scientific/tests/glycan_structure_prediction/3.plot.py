#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  simple_cycpep_predict/3.plot.py
## @brief This script is part of simple_cycpep_predict scientific test
## @author Sergey Lyskov.
## @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

import os, sys, subprocess, math, copy
import matplotlib
#matplotlib.use('Agg')

import matplotlib as mpl
import matplotlib.pyplot as plt
import benchmark
import pandas, seaborn
import matplotlib.pylab as plb
from typing import *


bm_type="default"
mpl.rcParams.update(mpl.rcParamsDefault)

#=================================
def set_rc_params():
	mpl.rcParams['axes.xmargin'] = 0  # x margin.  See `axes.Axes.margins`
	mpl.rcParams['axes.ymargin'] = 0  # y margin See `axes.Axes.margins`

	mpl.rcParams['polaraxes.grid'] = True    # display grid on polar axes
	#axes3d.grid         : True    # display grid on 3d axes

	### TICKS
	# see http://matplotlib.org/api/axis_api.html#matplotlib.axis.Tick
	mpl.rcParams['xtick.major.size'] = 4      # major tick size in points
	mpl.rcParams['xtick.minor.size'] = 2      # minor tick size in points
	mpl.rcParams['xtick.major.width'] = 1    # major tick width in points
	mpl.rcParams['xtick.minor.width'] = 1    # minor tick width in points
	mpl.rcParams['xtick.major.pad'] = 6      # distance to major tick label in points
	mpl.rcParams['xtick.direction'] = 'out'     # direction: in, out, or inout

	mpl.rcParams['ytick.major.size'] = 4      # major tick size in points
	mpl.rcParams['ytick.minor.size'] = 2      # minor tick size in points
	mpl.rcParams['ytick.major.width'] = 1    # major tick width in points
	mpl.rcParams['ytick.minor.width'] = 1    # minor tick width in points
	mpl.rcParams['ytick.major.pad'] = 6      # distance to major tick label in points
	#ytick.minor.pad      : 4      # distance to the minor tick label in points
	#ytick.color          : k      # color of the tick labels
	#ytick.labelsize      : medium # fontsize of the tick labels
	mpl.rcParams['ytick.direction'] = 'out'     # direction: in, out, or inout
	mpl.rcParams['figure.figsize'] = ( 7.5, 6)    # figure size in inches]

#=================================
def set_subplot_ylim_offset(gg, pad_percent=.05):
	for ax in gg.axes:
		#print(ax.get_ylim()[0], ax.get_ylim()[1])
		current = ax.get_ylim()[0]
		#current_upper = ax.get_ylim()[]
		d = abs(current - ax.get_ylim()[1])
		p = d * pad_percent
		#print("D:",  d)
		#print("P: ", p)

		new_lim = current - p
		new_lim_upper = ax.get_ylim()[1] + p
		#print("Old", current)
		#print("New", new_lim)
		#print("NewU", new_lim_upper)
		ax.set_ylim(new_lim, new_lim_upper)


#=================================
def plot_score_vs_rmsds(local_df: pandas.DataFrame, plots: List[Tuple[str,str]]):
	"""
	Plot score vs RMSD plots of top percent of models
	"""
	top = [.8, .15]

	if debug:
		top = [1]

	for top_p in top:

		top_n = (len(local_df) / len(local_df[['experiment', 'pdb_branch']].drop_duplicates())) * top_p
		print(top_n)

		clean = local_df.sort_values(['experiment', 'pdb_branch', 'total_score']).groupby(['experiment', 'pdb_branch']).head(
			top_n)
		print(len(clean))

		for name, f in clean.groupby(['experiment'], as_index=False):
			print(name)

			# if (name != "rounds_1"): continue
			print(name, len(f))
			if len(f) == 0: continue

			print("fit.6")
			g = seaborn.FacetGrid(f, col="pdb-branch-low-size", col_wrap=5, sharey=False, sharex=False)
			g.map(plb.scatter, "fit6_rmsd", "total_score", s=3 ** 2)
			g.set_titles("{col_name}")
			g.set_axis_labels("RMSD (fit.6)", "Total Score (Top " + str(top_p * 100) + "%)")
			g.set(xlim=(0, plt.xlim()[1]+1))
			set_subplot_ylim_offset(g, .1)

			title = "Score vs RMSD, Top "+str(top_p*100)+"% "+name
			out = bm_type + "_cleaned_" + str(top_p) + "_rmsd_vs_score_fit6_" + name + ".png"

			title_path = (title, out)
			plots.append(title_path)
			g.savefig(out, dpi=300)
			# plt.show()
			plt.close()

#=================================
def plot_kdes(local_df: pandas.DataFrame, plots: List[Tuple[str,str]]):
	"""
	Plot KDE of RMSD at different top energies
	"""

	top = [1, .25, .05, .01]
	if debug:
		top = [1, .75, .25]

	for top_p in top:
		top_n = (len(local_df) / len(local_df[['experiment', 'pdb_branch']].drop_duplicates())) * top_p
		print(top_n)

		clean = local_df.sort_values(['experiment', 'pdb_branch', 'total_score']).groupby(
			['experiment', 'pdb_branch']).head(top_n)
		print(len(clean))

		clean.groupby('experiment')['fit6_rmsd'].plot(kind='kde', x='fit6_rmsd', label="Top "+str(top_p*100)+"%")


		plt.xlabel('RMSD (fit.6)')

		plt.xlim(0, 30)
		plt.ylim(0, plt.ylim()[1] + .01)
		plt.ylim()
		out = "rmsd_density_distributions.png"

	title = "RMSD Distribution"
	plt.title(title)
	plt.legend()
#	title_path = (title, out)
#	plots.append(title_path)
#	plt.savefig(out, dpi=300)
		#plt.show()
#	plt.close()

#=================================
def plot_lowest_energy_rmsd(local_df: pandas.DataFrame, plots: List[Tuple[str,str]]):
	"""
	Plot the (fit6) RMSD of the lowest energy model for each input structure.
	:param local_df:
	:param plots:
	:return:
	"""
	top_scoring = local_df.reset_index().sort_values(['experiment', 'pdb_branch', 'total_score']).groupby(
		['experiment', 'pdb_branch']).head(1)
	top_scoring = top_scoring.sort_values('rmsd')
	top_scoring = top_scoring.sort_values('fit6_rmsd')
	seaborn.barplot(data=top_scoring, x='pdb-branch-low-size', y='fit6_rmsd', hue='experiment', hue_order=exp_order)
	title = "Lowest Energy Model \n >.6 Residue Fit"
	plt.ylim(0, 30)
	plt.title( title, fontsize=12)
	plt.ylabel("Heavy Atom RMSD", fontsize=10)
	plt.tick_params('x', labelrotation=90, labelsize=8)
	plt.xlabel("")
	plt.legend(fontsize=8)
	plt.gcf().subplots_adjust(bottom=0.25)
	out = bm_type + "_per_glycan_rmsd_fit6.png"

	title_path = (title, out)
	plots.append(title_path)
	plt.savefig(out, dpi=300)
	#plt.show()
#	plt.close()

#=================================
benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

## Load Dataframes ##
df = pandas.read_csv("data.csv")
enrichment_df = pandas.read_csv("sampling_enrichments.csv")

## Setup Matplotlib  ##
set_rc_params()

#width = 7.5
#height = 6
#plt.rc("font", size=12)
#plt.rcParams['figure.figsize'] = width, height #width, height


exp_order = sorted(df['experiment'].drop_duplicates())
plot_paths = []


## Do the plotting ##
plot_score_vs_rmsds(df, plot_paths)

mpl.rcParams['figure.figsize'] = ( 8, 5)    # figure size in inches

plt.subplot( 1, 2, 1 )
plot_kdes(df, plot_paths)

plt.subplot( 1, 2, 2 )
plot_lowest_energy_rmsd(df, plot_paths)


benchmark.save_variables('debug working_dir testname plot_paths failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
