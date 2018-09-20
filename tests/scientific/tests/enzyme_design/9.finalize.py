#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  enzyme_design/9.finalize.py
## @brief this script is part of enzyme_design scientific test
## @author Rocco Moretti



import json
import benchmark
import operator

from benchmark import *
from benchmark.tests import _PlotsKey_, _TestsKey_

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals

_index_html_template_ = '''\
<html>
<head>
    <title>Enzyme Design scientific benchmark test results</title>
  <style>
    fixed {{background-color: #eee; white-space: pre-wrap; font-family: Monaco, 'Liberation Mono', Courier, monospace; font-size:12px; }}
  </style>
</head>
<body>
<H2>Enzyme Design Benchmark</H2>
<p>
    The enzdes benchmark is described in <a href="https://doi.org/10.1002/prot.24463">Nivon et al. (2014)</a> "Automating human intuition for protein design."
    The concept is to redesign the binding sites of a set of proteins which bind to their native ligands.
    As the sequences of the native proteins should theoretically be optimal for binding, an ideal design protocol should optimize sequence recovery.

    The thresholds for this test are somewhat arbitrary, and are based off of the previous runs for this work.
</p>
<br/>
<h3>FAILURES</h3>
<p>
    {failures}
</p>
<h3>RESULTS</h3>
<p>
    {results}
</p>
<br/>
</body></html>
'''

def check_value( value, comparison, threshold ):
    funcs = {
        "<": operator.lt,
        ">": operator.gt,
        "<=": operator.le,
        ">=": operator.ge,
        "!=": operator.ne,
        "=": operator.eq,
        "==": operator.eq,
    }
    return funcs[comparison](value,threshold)

def failures(results, thresholds):
    failures = [ ]
    for target in thresholds:
        for metric, (comparison, value) in thresholds[target].items():
            if target not in results or metric not in results[target]:
                failures.append( ( target, metric, "Needs to be present" ) )
                continue
            result = results[target][metric]
            if not check_value(result, comparison, value):
                failures.append( ( target, metric, f"Wanted {result:.3f} {comparison} {value}" ) )
    return failures

def format_failures( fails ):
    output = ["<ul>"]
    for (target, metric, reason) in fails:
        output.append( "<li>{target} with metric {metric}: {reason}</li>" )
    if len(fails) == 0:
        output.append( "<i>(NONE)</i>" )
    output.append("</ul>")
    return '\n'.join(output)

def format_results_row( target, row, metrics ):
    output = []
    output.append( "<tr>" )
    output.append( f"<td>{target}</td>" )
    for metric in metrics:
        if metric in row:
            value = "{:.3f}".format( row[metric] )
        else:
            value = "N/C"
        output.append( f'<td align="center">{value}</td>' )
    output.append( "</tr>" )

    return '\n'.join( output )

def format_results( results ):
    output = ['<table>']
    metrics = results['ALL'].keys()

    output.append( "<tr>" )
    output.append( "<th>Protein</th>" )
    for metric in metrics:
        output.append( f'<th align="center">{metric}</th>' )
    output.append( "</tr>" )

    output.append( format_results_row( "AVERAGE", results['ALL'], metrics ) )

    targets = [ k for k in results.keys() if k != "ALL" ]
    targets.sort()
    for target in targets:
        output.append( format_results_row( target, results[target], metrics ) )

    output.append( "</table>" )

    return '\n'.join(output)

def write_html(results,thresholds) :
    fails = [ (f[0]!="ALL", f) for f in failures(results,thresholds) ]
    fails.sort()
    fails = [ b for (a,b) in fails ]

    with open(f'{working_dir}/index.html', 'w') as f:
        f.write(
            _index_html_template_.format(
                results=format_results( results ),
                failures=format_failures( fails ),
            )
        )

    if len(fails) > 0:
        return _S_failed_
    else:
        return _S_passed_

with open(_multi_step_result_, 'w') as f:
    overall_result = write_html(results, thresholds)
    r = {
        _StateKey_  : overall_result,
        _ResultsKey_ : {
            _PlotsKey_ : [
                {
                    "data": [
                      {
                        "color": "#66f",
                        "legend": "Sequence Recovery",
                        "y": "seqrec"
                      }
                    ],
                    "type": "sub_test:revision_value"
                },
                {
                    "data": [
                      {
                        "color": "#66f",
                        "legend": "Fraction PSSM Recovery",
                        "y": "pssm_seqrec"
                      }
                    ],
                    "type": "sub_test:revision_value"
                },
                {
                    "data": [
                      {
                        "color": "#66f",
                        "legend": "Average PSSM Score Delta",
                        "y": "pssm_delta_seqrec"
                      }
                    ],
                    "type": "sub_test:revision_value"
                },
            ],
            _TestsKey_: {
                "OVERALL_AVERAGE": {
                    "seqrec": results["ALL"]["seqrec"],
                    "pssm_seqrec": results["ALL"]["pssm_seqrec"],
                    "pssm_delta_seqrec": results["ALL"]["pssm_delta_seqrec"],
                    _StateKey_: overall_result,
                    _LogKey_ : "Average results from {} structures.".format(len(targets))
                },
            },
        },
        _LogKey_ : 'Done!',

    }

    json.dump(r, f, sort_keys=True, indent=2)

benchmark.save_variables()  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
