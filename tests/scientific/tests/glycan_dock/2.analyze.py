#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  glycan_dock/2.analyze.py
## @brief This script is part of the GlycanDock scientific test
## @author Sergey Lyskov
## @author Morgan Nance (@mlnance; revised for GlycanDock sci test)

import os, sys, subprocess, math
import numpy as np
import benchmark
from benchmark.util import quality_measures as qm
from benchmark.util import bootstrap as bs
from benchmark.util import scorefile_io_utils as sciu


benchmark.load_variables()  # Python black magic: load all variables saved by previous script into s
config = benchmark.config()

results = {}
scorefiles = []
#logfiles = []
cutoffs_ring_Lrmsd_dict = {}
cutoffs_bootstrap_N5_dict = {}
failures = []

# inputs are header labels from the scorefile
# example: "total_score" and "rmsd"
x_label = "ring_Lrmsd"
x_cutoff = 2.0 # near-native model is <= 2.0 ring_Lrmsd
y_label = "interaction_energy"
outfile = "result.txt"
cutoffs = "cutoffs"

# scorefiles and logfiles
scorefiles.extend( [ f'{working_dir}/output/{t}/{t}.score' for t in targets ] )

# read cutoffs
protein = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_ring_Lrmsd = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $2}'" ).splitlines()
cutoffs_bootstrap_N5 = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $3}'" ).splitlines()
cutoffs_ring_Lrmsd = map( float, cutoffs_ring_Lrmsd )
cutoffs_bootstrap_N5 = map( float, cutoffs_bootstrap_N5 )
cutoffs_ring_Lrmsd_dict.update( dict( zip ( protein, cutoffs_ring_Lrmsd )))
cutoffs_bootstrap_N5_dict.update( dict( zip ( protein, cutoffs_bootstrap_N5 )))

# open results output file
fh = open( outfile, "w" )

# go through scorefiles of targets
topX = 5
for sc_f, target in zip(scorefiles, targets):

        target_results = {}

        # read in scorefile as a pandas DataFrame
        df = sciu.scorefile_to_dataframe(sc_f)
        df.to_csv(os.path.join(os.path.dirname(sc_f), "{}.csv".format(target)), index=None)
        # sort data on y_label i.e. interaction_energy
        # and renumber so that index goes from 0, 1, 2, etc
        df = df.sort_values(y_label).reset_index(drop=True)

        # calculate the bootstrap average N5 sorted by interaction_energy
        # look at the top-5 scoring by y_label and see how many are
        # docking successes, where success is x_label <= x_cutoff
        # False here means gt_eq_cutoff = False i.e. <= x_cutoff
        nscore, nstd = bs.bootstrap_NX( sc_f, y_label, 5, x_label, x_cutoff, False)

        # check for bootstrap_N5 >= cutoff
        # unique bootstrap_N5 (<N5>) cutoff per target
        fh.write( target + "\t" )
        bs_tag = "<N5>"
        val_cutoff = qm.check_avgNX_above_cutoff(nscore, cutoffs_bootstrap_N5_dict[target], bs_tag, fh)
        target_results.update( val_cutoff )
        # ONLY FAILURE case if bootstrap_N5 not >= cutoff
        if val_cutoff["Average " + bs_tag + " >= cutoff"] == False and target not in failures:
            failures.append( target )

        # report <N5> range
        fh.write( target + "\t" )
        val_N5 = qm.get_N5_spread(nscore, nstd, bs_tag, fh)
        target_results.update( val_N5 )

        # report if 1 or more of top-5 models is <= its cutoff ring_Lrmsd
        fh.write( target + "\t" )
        val_topXscoring = qm.check_rmsd_of_topXscoring( df[x_label], topX, cutoffs_ring_Lrmsd_dict[target], fh )
        target_results.update( val_topXscoring )

        # report full ring_Lrmsd range
        fh.write( target + "\t" )
        val_rmsd = qm.check_range( df[x_label], x_label, fh )
        target_results.update( val_rmsd )

        # report 5-top-scoring ring_Lrmsd range
        fh.write( target + "\t" )
        val_rmsd_topX = qm.check_range( df[x_label].iloc[:topX], x_label + " of top-" + str(topX), fh )
        target_results.update( val_rmsd_topX )

        # report full interaction_energy range
        fh.write( target + "\t" )
        val_score = qm.check_range( df[y_label], y_label, fh )
        target_results.update( val_score )

        # report 5-top-scoring interaction_energy range
        fh.write( target + "\t" )
        val_score_topX = qm.check_range( df[y_label].iloc[:topX], y_label + " of top-" + str(topX), fh )
        target_results.update( val_score_topX )

        # and report PNear just for info
        pnear = qm.calculate_pnear( df[y_label], df[x_label], lambda_val=2.0, kbt=0.62 )
        target_results.update( {"PNear (lambda = 2.0)": round(pnear, 3)} )

        results.update( {target : target_results} )
        fh.write( "\n" )

fh.close()

benchmark.save_variables('debug targets nstruct working_dir testname results scorefiles cutoffs_ring_Lrmsd_dict cutoffs_bootstrap_N5_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
