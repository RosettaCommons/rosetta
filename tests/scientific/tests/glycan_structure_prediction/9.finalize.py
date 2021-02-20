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

overall_pass = True
if len(failures) > 0:
    overall_pass = False


# build up html from readme, start with the starting tag
_index_html_template_ = "<html>\n"
_index_html_template_ += "<H2>Scientific test: " + testname + "</H2>\n"
_index_html_template_ += "<h3>FAILURES</h3>\n<p>\n"

# add failures to html
if len( failures ) > 0:
    for failure in failures:
        _index_html_template_ += str(failure) + "<br>\n"

_index_html_template_ += "<h3>RESULTS</h3>\n"


# add image for results
for title_path in plot_paths:
    title = title_path[0]
    path = title_path[1]
    _index_html_template_ += f'<img src="{path}" alt="{title}" style="max-width: 100%">\n'

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
    if overall_pass == True or debug:
        overall_result = _S_passed_
    else:
        overall_result = _S_failed_
    write_html( )
    r = {
        _StateKey_  : overall_result,
        _ResultsKey_ : {},
        _LogKey_ : 'Done!',
    }

    json.dump(r, f, sort_keys=True, indent=2)

benchmark.save_variables()  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
