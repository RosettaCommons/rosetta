#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  stepwise_rna_favorites/9.finalize.py
## @brief this script is part of stepwise_rna_favorites scientific test
## @author Andy Watkins


import json
import benchmark

from benchmark import *

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals

def failures(results):
    s = "<b>Failures:</b><br/><ol>"
    at_least_one = False
    for target in targets:
        if "FALSE" in list(results[target].values()):
            s += "<li>" + target + "</li>"
            at_least_one = True
    s += "</ol>"
    if at_least_one: return s
    else: return "<b>No failures; test passed.</b>"

_index_html_template_ = '''\
<html>
<head>
    <title>SWM scientific benchamrk test results</title>

  <style>
    fixed {{background-color: #eee; white-space: pre-wrap; font-family: Monaco, 'Liberation Mono', Courier, monospace; font-size:12px; }}
  </style>
</head>
<body>
<p>
    The SWM benchmark, run on the classic RNA benchmark set 'favorites.txt', evaluates the stability in scientific performance of
    the Das lab's high resolution Stepwise Monte Carlo structure prediction protocol. Specifically, we run twelve 'motif-scale' structure prediction targets of varying difficulty
    and gauge whether the minimum energy structure's RMSD and the maximum sampled RMSD are within appropriate values.
</p>
<br/>
<p>
    {failures}
</p>
<br/>
<p>
    <img src='score_vs_rmsd.png' width='1200' height='720' />
</p>
</body></html>
'''


def state_of_results(results):
    """
    Judges whether the test has passed given this dict of targets-to-stats.
    """
    with open(f'{working_dir}/index.html', 'w') as f:

        f.write(
            _index_html_template_.format(
                failures=failures(results),
            )
        )

    for key, value in results.items():
        if "FALSE" in value.values():
            return _S_failed_

    return _S_passed_

with open(_multi_step_result_, 'w') as f:
    r = {
        _StateKey_  : state_of_results(results),
        _ResultsKey_ : {},
        _LogKey_ : 'Done!',
    }

    json.dump(r, f, sort_keys=True, indent=2)


benchmark.save_variables()  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
