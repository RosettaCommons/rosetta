#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  ligand_docking/9.finalize.py
## @brief this script is part of the ligand docking scientific test
## @author Sergey Lyskov


import json, subprocess
import benchmark

from benchmark import *

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

cutoff = 1.0 ## Due to change once we run things a couple of times.
# read readme
with open('readme.md') as f: readme = f.read().splitlines()

# build up html from readme, start with the starting tag
_index_html_template_ = "<html>\n"
_index_html_template_ += "<H2>Scientific test: " + testname + "</H2>\n"
_index_html_template_ += "<h3>SAMPLING FAILURES</h3>\n<p>\n"

if len( sampling_failures ) > 0:
    for sf in sampling_failures.keys():
        _index_html_template_ += str(sf) + ": " + str(len(sampling_failures[sf])) + "<br>\n"
        for target in sampling_failures[sf]:
            _index_html_template_ += str(target) + "\t"
        _index_html_template_ += "<br>\n"

else:
    _index_html_template_ += "&nbsp;&nbsp;&nbsp;&nbsp;None<br>\n"

##scoring failures
_index_html_template_ += "<h3>SCORING FAILURES</h3>\n<p>\n"

if len( scoring_failures ) > 0:
    for sf in scoring_failures.keys():
        _index_html_template_ += str(sf) + ": " + str(len(scoring_failures[sf])) + "<br>\n"
        for target in scoring_failures[sf]:
            _index_html_template_ += str(target) + "\t"
        _index_html_template_ += "<br>\n"
_index_html_template_ += "</p>\n<h3>RESULTS</h3>\n"
_index_html_template_ += '<img src="plot_results.png" alt="alternative text" style="max-width: 100%">\n'

# add text from readme
for l in readme:

    # headings
    if l.startswith( "## " ):
        _index_html_template_ += "<h3>" + l.replace( ">> ", "" ) + "</h3>\n"
    
    # ignore the description
    elif l.startswith( "#### " ):
        continue
        
    # insert the actual text as a paragraph
    else:
        _index_html_template_ += "<p>" + l + "</p>\n"

_index_html_template_ += "</body></html>\n"

## Looking at worst performing score function to compare to cutoff value.
scoring_failures_by_sfxn = []
for sfxn in sfxns:
    scoring_failures_by_sfxn.append(len(scoring_failures[sfxn]))
for_cutoff = float(max(scoring_failures_by_sfxn)/len(targets))

# write the html
def write_html( sampling_failures , scoring_failures ) :
    with open(f'{working_dir}/index.html', 'w') as f:
        f.write( _index_html_template_.format( sampling_failures ) )

    if for_cutoff <= cutoff or config['debug']:
        return _S_passed_
    else:
        return _S_failed_

# write the overall results
with open(_multi_step_result_, 'w') as f:
    overall_result = write_html( sampling_failures, scoring_failures )
    r = {
        _StateKey_  : overall_result,
        _ResultsKey_ : {},
        _LogKey_ : 'Done!',
    }

    json.dump(r, f, sort_keys=True, indent=2)


benchmark.save_variables()  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
