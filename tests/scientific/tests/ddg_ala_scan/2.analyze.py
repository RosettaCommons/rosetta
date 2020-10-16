#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/2.analyze.py
## @brief this script is part of cartesian_relax scientific test
## @author Sergey Lyskov

import os, sys, subprocess, math

import json

import pandas as pd
import sqlite3
import itertools
from glob import glob
import benchmark
from benchmark.util import quality_measures as qm
import numpy as np

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	
#debug targets nstruct working_dir testname
config = benchmark.config()

#============================================
def parse_db_output(db_output_file, ddg_data):
    """Function to parse ddg calculated data from dbe sqlite files. Files created by ReportToDB Mover"""
    query = 'SELECT ddG.resNum, ddg.mutated_to_name3, ddg.ddG_value, ' \
    'structures.tag, residue_pdb_identification.residue_number, ' \
    'residue_pdb_identification.chain_id, residues.name3, batches.description, batches.name ' \
    'FROM ddg INNER JOIN structures ON structures.struct_id=ddg.struct_id ' \
    'INNER JOIN residue_pdb_identification ON residue_pdb_identification.struct_id=structures.struct_id ' \
    'AND residue_pdb_identification.residue_number=ddg.resNum ' \
    'INNER JOIN residues ON residues.struct_id=structures.struct_id AND residues.resNum=ddg.resNum ' \
    'INNER JOIN batches ON batches.batch_id=structures.batch_id'

    conn = sqlite3.connect(db_output_file)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    for rosetta_resnum, mutated_to, ddg_calc, tag, pdb_resnum, chain, original_name3, description, PDBPosID in c.execute(query):
        pdb_id = tag.split('_')[0]
        row = {}
        row['PDBPosID']=PDBPosID 
        row['score_function'] = description
        row['DDG_calc'] = ddg_calc 
        ddg_data.append(row)

    conn.close()
    return ddg_data

#============================================
def get_fraction_correctly_classified(df, e_cutoff = 1):
    """Returns the fraction of mutations that are correctl classified. 
    So both experimental and calculated ddg are either neutral, both positive or both negative.
    Expects a pandas dtaframe with columns DDG_obs and DDG_calc"""

    num_both_positive = len(ddg_df.query(f'DDG_obs>{e_cutoff} and DDG_calc>{e_cutoff}'))
    num_both_negative = len(ddg_df.query(f'DDG_obs<{-e_cutoff} and DDG_calc<{-e_cutoff}'))
    num_both_neutral  = len(ddg_df.query(f'abs(DDG_obs)<{e_cutoff} and abs(DDG_calc)<{e_cutoff}'))
    all_mut = len(ddg_df)
    fraction_correct_classified = (num_both_positive + num_both_negative + num_both_neutral)/ all_mut
    return fraction_correct_classified

#============================================

#Load calculated ddgs
db_files = glob(f"{working_dir}/output/*/*.db3")
ddg_data = []
for db_file in db_files:
    parse_db_output(db_file,ddg_data)

ddg_calc_df = pd.DataFrame(ddg_data)
ddg_calc_df.set_index('PDBPosID', inplace=True)

#Load experimental data
ddg_exp_df = pd.read_csv('mutation_benchmark_set.csv')
#some entries say >xxx, i.e. specify only a lower bond. Do we want to use this points?
ddg_exp_df['lower_bound']=ddg_exp_df.DDG_obs.str.contains('>')
ddg_exp_df.DDG_obs=ddg_exp_df.DDG_obs.str.replace('>','')
ddg_exp_df.DDG_obs=ddg_exp_df.DDG_obs.astype(float)
ddg_exp_df.drop(columns='DDG_calc', inplace=True)
ddg_exp_df.set_index('PDBPosID', inplace=True, drop=False)


results = {}
failures = []

#Join experimental and calculated data on PDBPosID
ddg_df = ddg_exp_df.join(ddg_calc_df, how='left' )

ddg_df.to_csv(f"{working_dir}/output/calc_vs_exp.csv")


score_function = 'talaris2014'

#fit the correlation coefficient
from scipy.stats import linregress
import numpy as np

fit = linregress(ddg_df.DDG_obs,ddg_df.DDG_calc)
MAE = np.mean(np.abs(ddg_df.DDG_obs-ddg_df.DDG_calc))
FCC = get_fraction_correctly_classified(ddg_df)

print("Slope: ", fit.slope)
print("intercept: ", fit.intercept)
print("rvalue: ", fit.rvalue)
print("pvalue: ", fit.pvalue)
print("stderr: ", fit.stderr)

# read cutoffs and put that into a dictionary
measures = subprocess.getoutput( "awk '{print $1}' cutoffs" ).splitlines()
cutoffs = subprocess.getoutput( "awk '{print $2}' cutoffs" ).splitlines()
cutoff_dict = dict( zip( measures, cutoffs ) )

results[score_function] = {'fit':fit, 'R':fit.rvalue, 'MAE':MAE, 'FCC':FCC}

# compare and get failures
min_R = float( cutoff_dict['min_R'] )
if fit.rvalue < min_R:
    failures += [f'{score_function}: R is too low! Should be higher than {min_R} but is {fit.rvalue}.']

max_MAE = float( cutoff_dict['max_MAE'] )
if MAE > max_MAE:
    failures += [f'{score_function}: MAE (MEAN ABSOLUTE ERROR) is too High! Should be lower than {max_MAE} but is {MAE}.']

min_FCC = float( cutoff_dict['min_FCC'] )
if FCC < min_FCC:
    failures += [f'{score_function}: FCC (Fraction correctly classified) is too low! Should be higher than {min_FCC} but is {FCC}.']

benchmark.save_variables('debug targets nstruct working_dir testname results  failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
