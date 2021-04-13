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

cutoff = {'ligand': 33, 'talaris':28, 'ref2015':29, 'betanov16':34}

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
_index_html_template_ += '<img src="plot_results1.png" alt="alternative text" style="max-width: 100%">\n'
_index_html_template_ += '<img src="plot_results2.png" alt="alternative text" style="max-width: 100%">\n'
_index_html_template_ += '<img src="plot_results3.png" alt="alternative text" style="max-width: 100%">\n'
_index_html_template_ += '<img src="plot_results4.png" alt="alternative text" style="max-width: 100%">\n'
_index_html_template_ += '<img src="plot_results5.png" alt="alternative text" style="max-width: 100%">\n'
_index_html_template_ += '<img src="plot_results6.png" alt="alternative text" style="max-width: 100%">\n'
_index_html_template_ += '<img src="plot_results7.png" alt="alternative text" style="max-width: 100%">\n'

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
number_of_sfxn_failures = 0
for sfxn in sfxns:
	if scoring_failures[sfxn] > cutoff[sfxn]:
		number_of_sfxn_failures += 1

# write the html
def write_html( sampling_failures , scoring_failures ) :
    with open(f'{working_dir}/index.html', 'w') as f:
        f.write( _index_html_template_.format( sampling_failures ) )

    if number_of_sfxn_failures == 0 or config['debug']:
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
