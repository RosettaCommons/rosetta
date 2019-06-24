#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mhc_epitope_energy/9.finalize.py
## @brief this script is part of mhc_epitope_energy scientific test
## @author Sergey Lyskov
## @author Brahm Yachnin (brahm.yachnin@rutgers.edu)


import json, subprocess
import benchmark

from benchmark import *
from benchmark.tests import _TestsKey_

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

# read readme
readme = subprocess.getoutput( "cat readme.md" ).splitlines()

# build up html from readme, start with the starting tag
_index_html_template_ = "<html>\n"
_index_html_template_ += "<H2>Scientific test: " + testname + "</H2>\n"
_index_html_template_ += "<h3>FAILURES</h3>\n<p>\n"

# Dictionary of subtests, keeping track of the results.  Assume everything passed.
subtest_dict = {key:{_StateKey_:_S_passed_, _LogKey_:"Successfully ran " + key + " test with " + str(len(targets)) + " targets."} for key in ["mhc_epitope", "delta_mhc_epitope", "base_total_score_vs_mhc_epitope", "sequence_recovery_vs_mhc_epitope", "core_sequence_recovery_vs_mhc_epitope", "delta_packstat_vs_mhc_epitope", "delta_buried_unsat_vs_mhc_epitope", "delta_netcharge_vs_mhc_epitope"]}

# add failures to html, and keep track of failed subtests in the subtest_dict
if len( failures ) > 0:
    if config['debug']: _index_html_template_ += "mhc_epitope_energy was run in DEBUG mode.  These test failures are not meaningful.<br>\n"
    for failedpdb, failedtests in failures_dict.items():
        _index_html_template_ += str(failedpdb) + ": " + ", ".join(failedtests) + "<br>\n"
        for subtest in failedtests:
            subtest_dict.update({subtest:{_StateKey_:_S_failed_, _LogKey_:subtest + " test ran with failures!"}})
else:
    _index_html_template_ += "&nbsp;&nbsp;&nbsp;&nbsp;None<br>\n"
	
_index_html_template_ += "</p>\n<h3>README</h3>\n"

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
		
_index_html_template_ += "</p>\n<h3>RESULTS</h3>\n"

#List of results
for result in subtest_dict.keys():
    _index_html_template_ += "</p>\n<h4>" + result + "</h4>\n"
    _index_html_template_ += '<img src="plot_results_' + result + '.png" alt="alternative text" style="max-width: 100%">\n'

# html closing tag
_index_html_template_ += "</body></html>\n"


# write the html
def write_html( failures ) :
    with open(f'{working_dir}/index.html', 'w') as f:
        f.write( _index_html_template_.format( failures ) )

    if len(failures) == 0 or config['debug']:
        return _S_passed_
    else:
        return _S_failed_

# write the overall results
with open(_multi_step_result_, 'w') as f:
    overall_result = write_html( failures )
    r = {
        _StateKey_  : overall_result,
        _ResultsKey_ : {
			_TestsKey_ : subtest_dict
		},
        _LogKey_ : 'Done!',
    }

    json.dump(r, f, sort_keys=True, indent=2)


benchmark.save_variables()  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
