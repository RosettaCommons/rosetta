#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  simple_cycpep_predict/9.finalize.py
## @brief this script is part of simple_cycpep_predict scientific test
## @author Sergey Lyskov
## @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).


import json, subprocess
import benchmark

from benchmark import *

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

# read readme
readme = subprocess.getoutput( "cat readme.md" ).splitlines()

# build up html from readme, start with the starting tag
_index_html_template_ = "<html>\n"
_index_html_template_ += "<H2>Scientific test: " + testname + "</H2>\n"
_index_html_template_ += "<h3>FAILURES</h3>\n<p>\n"
if overall_pass == False:
    _index_html_template_ += "See results summary, below, for failures." + "<br/>\n"
else:
    _index_html_template_ += "None." + "<br/>\n"
_index_html_template_ += "<h3>RESULTS SUMMARY</h3>\n<p>\n"

# add pass/fail summary to html
summary_text = subprocess.getoutput( f"cat {working_dir}/result.txt" ).splitlines()
for line in summary_text:
    _index_html_template_ += line + "<br/>\n"

# add image for results
_index_html_template_ += '<img src="plot_results.png" alt="Energy landscape computed for peptide to known crystal structure.  This should look funnel-shaped." style="max-width: 100%">\n'
_index_html_template_ += '<img src="plot_results2.png" alt="Energy landscape computed for peptide to lowest-energy structure.  This should also look funnel-shaped." style="max-width: 100%">\n'

# add text from readme
for l in readme:

    # headings
    if l.startswith( "## " ):
        _index_html_template_ += "<h3>" + l.replace( "## ", "" ) + "</h3>\n"
    
    # ignore the description
    elif l.startswith( "#### " ):
        _index_html_template_ += "<p><i>" + l.replace( "#### ", "" ) + "</i></p>\n"
        
    # insert the actual text as a paragraph
    else:
        _index_html_template_ += "<p>" + l + "</p>\n"

# html closing tag
_index_html_template_ += "</body></html>\n"


# write the html
def write_html( ) :
    with open(f'{working_dir}/index.html', 'w') as f:
        f.write( _index_html_template_.format( ) )

# write the overall results
with open(_multi_step_result_, 'w') as f:
    if overall_pass == False:
        overall_result = _S_failed_
    else:
        overall_result = _S_passed_
    write_html( )
    r = {
        _StateKey_  : overall_result,
        _ResultsKey_ : {},
        _LogKey_ : 'Done!',
    }

    json.dump(r, f, sort_keys=True, indent=2)

benchmark.save_variables()  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
