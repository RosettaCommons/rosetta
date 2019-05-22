#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cofactor_binding_sites/9.finalize.py
## @brief this script is part of cofactor_binding_sites scientific test
## @author Amanda Loshbaugh (aloshbau@gmail.com)

import json, subprocess
import benchmark

from benchmark import *

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

r = {
    _ResultsKey_ : results,
    _LogKey_ : 'Done!',
}
if failure == False:
    r[_StateKey_] = _S_passed_
    r[_LogKey_] = 'The test passed! You can learn more about CoupledMoves at https://www.rosettacommons.org/docs/latest/coupled-moves'
elif failure == True:

    r[_StateKey_] = _S_failed_
    r[_LogKey_] = 'The test failed! Here are some things that might have gone wrong. (1) Did you change code for CoupledMoves, BoltzmannRotamerMover, or Backrub mover? CoupledMoves uses the BoltzmannRotamerMover for side chain moves, and the Backrub mover for backbone moves. (2) Did CoupledMoves run correctly? Analysis requires fasta files which are written by CoupledMoves when it runs. Check that those exist. For a production run, each fasta file should contain thousands of unique sequences. You can learn more about CoupledMoves at https://www.rosettacommons.org/docs/latest/coupled-moves'

# read readme
with open("readme.md") as f: readme = f.readlines()

# build up html from readme, start with the starting tag
_index_html_template_ = "<html>\n"
_index_html_template_ += "<H2>Scientific test: {0}</H2>\n".format( testname )
_index_html_template_ += "<h3>FAILURE:</h3>\n<p>\n"
_index_html_template_ += f"&nbsp;&nbsp;&nbsp;&nbsp;{failure}<br>\n"
_index_html_template_ += f"{r[_LogKey_]}<br>\n"#.format(  )
_index_html_template_ += "</p>\n<h3>RESULTS</h3>\n"
#_index_html_template_ += '<embed src="results.pdf" width="700" height="700">\n'
_index_html_template_ += '<img src="plot_results.png" alt="alternative text" style="max-width: 100%">\n'

# add text from readme
for l in readme:
    # headings
    if l.startswith( "## " ):
        _index_html_template_ += "<h3>{0}</h3>\n".format( l.replace( "## ", "" ) )
    # ignore the description
    elif l.startswith( "#### " ):
        continue
    # insert the actual text as a paragraph
    else:
        _index_html_template_ += "<p>{0}</p>\n".format( l )
# html closing tag
_index_html_template_ += "</body></html>\n"

# write the html
def write_html( failure, log ) :
    with open(f'{working_dir}/index.html', 'w') as f:
        f.write( _index_html_template_ )

f = open(_multi_step_result_, 'w')

json.dump(r, f, sort_keys=True, indent=2)

write_html( failure, r[_LogKey_] )

# write the overall results
# with open(_multi_step_result_, 'w') as f:
#     overall_result = write_html( failure )
#     r = {
#         _StateKey_  : overall_result,
#         _ResultsKey_ : {},
#         _LogKey_ : 'Done!',
#     }
#
#     json.dump(r, f, sort_keys=True, indent=2)


benchmark.save_variables()  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
