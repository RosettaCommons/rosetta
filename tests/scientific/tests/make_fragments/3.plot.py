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
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import benchmark

# benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()
debug = config['debug']


def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_ylabel('min rmsd fragment rmsd')

final_figure_file = "plot_results.png"

# TODO results file from script 2.
results_file = "2.checkpoint.json"
with open(results_file) as fh:
    results_data = json.loads(fh.read())

fig = plt.figure(constrained_layout=True)
spec = gridspec.GridSpec(ncols=10, nrows=10, figure=fig)
precision_plot = fig.add_subplot(spec[0:5, 0:7])
coverage_plot = fig.add_subplot(spec[5:, 0:7])
ax_scatter_plot = fig.add_subplot(spec[:, 7:])
# ax_summary = fig.add_subplot(spec[3:, 0])

precision_plot.set_xlabel("rmsd")
precision_plot.set_ylabel("Precision")

coverage_plot.set_xlabel("rmsd")
coverage_plot.set_ylabel("Coverage")

# TODO we can do this for everything, or plot every target at a time
# Right now this will only run once
datas = []
all_precisions = []
all_coverages = []
for rank_i, rank_limit in enumerate([10, 25, 50]):
    rmsd_cutoffs = np.linspace(0, 5, 60)
    precisions = []
    coverages = []

    for rmsd_cutoff in rmsd_cutoffs:
        totals = 0
        goods = 0
        coverage_count = 0
        coverage_total = 0
        for db_name, db_data in results_data.items():
            for target_name, target_data in db_data.items():
                for i, rmsd_set in enumerate(target_data['rmsds']):
                    if not rmsd_set:
                        continue
                    totals += rank_limit
                    count = len([x for x in rmsd_set[:rank_limit] if x <= rmsd_cutoff])
                    if min(rmsd_set[:rank_limit]) <= rmsd_cutoff:
                        coverage_count += 1
                    coverage_total += 1
                    goods += count
        precisions.append(goods/totals)
        coverages.append(coverage_count/coverage_total)
    all_precisions.append(precisions)
    all_coverages.append(coverages)

    mins = []
    for db_name, db_data in results_data.items():
        for target_name, target_data in db_data.items():
            for i, rmsd_set in enumerate(target_data['rmsds']):
                if not rmsd_set:
                    continue
                mins.append(min(rmsd_set[:rank_limit]))
    datas.append(mins)
    precision_plot.plot(rmsd_cutoffs, precisions, label=f"Precision {rank_limit}", marker='.')
    coverage_plot.plot(rmsd_cutoffs, coverages, label=f"Coverage {rank_limit}", marker='.')


parts = ax_scatter_plot.violinplot(datas, showmeans=False, showmedians=False, showextrema=False)

for pc in parts['bodies']:
    pc.set_facecolor('#D43F3A')
    pc.set_edgecolor('black')
    pc.set_alpha(1)

quartile1, medians, quartile3 = np.percentile(datas, [25, 50, 75], axis=1)

print(quartile1, medians, quartile3)
whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
            for sorted_array, q1, q3 in zip(datas, quartile1, quartile3)])
whiskersMin, whiskersMax = whiskers[:, 0], whiskers[:, 1]

inds = np.arange(1, len(medians) + 1)
ax_scatter_plot.scatter(inds, medians, marker='o', color='white', s=20, zorder=3)
ax_scatter_plot.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=4)
ax_scatter_plot.vlines(inds, whiskersMin, whiskersMax, color='k', linestyle='-', lw=1)

######### HERE we add results
failures = []
ret_dat = {
    "inds": inds,
    "quartile1": quartile1,
    "medians": medians,
    "quartile3": quartile3,
    "whistersMin": whiskersMin,
    "whistersMax": whiskersMax,
    "all_precisions": np.array(all_precisions),
    "all_coverages": np.array(all_coverages)}
ret_dat = {k: v.tolist() for k, v in ret_dat.items()}

with open("plot_check.json", 'w') as fh:
    fh.write(json.dumps(ret_dat))

expected_results_json_file = os.path.join(config["rosetta_dir"], "tests/scientific/data/make_fragments/expected_results.json")
with open(expected_results_json_file) as fh:
    all_expected_data = json.loads(fh.read())
if debug:
    expected_data = all_expected_data['debug']
else:
    expected_data = all_expected_data['full']

expected_medians = expected_data["medians"]

for i, result in enumerate(medians):
    if result > expected_medians[i] + 0.03:
        failures.append(f"Found that median of idx {i} was greater than 0.03 of the initial run: initial: {expected_medians[i]} current: {result}")

# set style for the axes
labels = ['<= Rank 10', '<= Rank 25', '<= Rank 50']
for ax in [ax_scatter_plot]:
    set_axis_style(ax, labels)

precision_plot.legend()
coverage_plot.legend()
fig.savefig(final_figure_file, dpi=300)

benchmark.save_variables('debug targets nstruct working_dir testname results outfile cutoffs_rmsd_dict cutoffs_score_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
