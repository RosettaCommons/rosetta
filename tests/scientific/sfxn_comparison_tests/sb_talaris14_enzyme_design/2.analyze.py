#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  enzyme_design/2.analyse.py
## @brief this script is part of enzyme_design scientific test
## @author Rocco Moretti

import os, sys, subprocess, math
import numpy as np
import benchmark
import json

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into
config = benchmark.config()

thresholds = {
    "ALL":{
        "seqrec":(">=",0.38),
    },
}

results = {}

results[ "ALL" ] = dict( (metric,0) for metric in metrics )

for target in targets: # Load targets from previous run.
    target_metric_list = dict( (metric,[]) for metric in metrics )
    with open( os.path.join( working_dir, "output", target, f"{target}.json" ) ) as f:
        for line in f: # The file format of the json-formatted scorefile is a per-line concatenated json
            data = json.loads(line)

            for metric in metrics:
                if metric in data:
                    target_metric_list[ metric ].append( data[metric] )

    target_metrics = {}
    for metric in metrics:
        if len( target_metric_list[ metric ] ) == 0:
            average = 0
        else:
            average = sum( target_metric_list[ metric ] ) / len( target_metric_list[ metric ] )

        target_metrics[ metric ] = average
        results["ALL"][ metric ] += average/len(targets)

    results[ target ] = target_metrics

benchmark.save_variables('targets nstruct working_dir testname results thresholds')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
