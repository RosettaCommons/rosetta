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
## @author Shannon Smith

import os, sys, subprocess, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import benchmark
import pandas as pd

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#=========================================
def plot_exp_vs_pred(data_dict):
    plt.xlabel("pEXP")
    plt.ylabel("Interface Score")

    for i, t in enumerate(data_dict.keys()):
        x_vals = []
        y_vals = []
        for l in data_dict.get(t):
            if data_dict[t][l]['interface_delta_X']['avg'] >= 10 :
                continue
            else:
                x_vals.append(float(data_dict[t][l]['log of experimental binding affinity ']))
                y_vals.append(float(data_dict[t][l]['interface_delta_X']['avg']))
        print(t, min(y_vals), max(y_vals))
        print(t, min(y_vals)*1.2, max(y_vals)/1.2)
        plt.xlim(int(math.floor(min(x_vals))), int(math.ceil(max(x_vals))))
        plt.ylim(min(y_vals)* 1.2 , int(max(y_vals) / 1.2))
        plt.plot(x_vals, y_vals, 'ko', markersize=4)
        plt.subplot(nrows, ncols, i + 1)
        plt.title(str(t) + ', SRC:' + str(spearman_rank_correlation(x_vals, y_vals)))

    ncols = 4
    nrows = 1
    if len(data_dict.keys()) < 4:
        ncols = len(data_dict.keys())
    else:
        nrows = math.ceil(len(data_dict.keys()) / 4)

    plt.subplots(ncols=ncols,nrows=nrows,figsize=(15,15))
    plt.subplots_adjust(wspace=0.2, hspace=0.75)
    plt.rc("font", size=8)

    plt.savefig('exp_vs_pred_ddG.png')


#=========================================

benchmark.save_variables('working_dir testname target_results outfile')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
plot_exp_vs_pred(results)
