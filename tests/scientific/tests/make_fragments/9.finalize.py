#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/9.finalize.py
## @brief this script is part of cartesian_relax scientific test
## @author Sergey Lyskov


import json, subprocess
import benchmark
import os

from benchmark import *

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()
debug       = config['debug']
testname = "make_fragments"
working_dir = os.getcwd()
print("working dir", working_dir)

# read readme
readme = subprocess.getoutput( "cat readme.md" ).splitlines()

# build up html from readme, start with the starting tag
_index_html_template_ = "<html>\n"
_index_html_template_ += "<H2>Scientific test: " + testname + "</H2>\n"
_index_html_template_ += "<h3>FAILURES</h3>\n<p>\n"

# add failures to html
if len( failures ) > 0:
    for failure in failures:
        _index_html_template_ += str(failure) + "<br>\n"
else:
    _index_html_template_ += "&nbsp;&nbsp;&nbsp;&nbsp;None<br>\n"
_index_html_template_ += "</p>\n<h3>RESULTS</h3>\n"
_index_html_template_ += '<img src="plot_results.png" alt="alternative text" style="max-width: 100%">\n'

in_list = False
list_lines = []

in_code = False
code_lines = []
# add text from readme
newline = "\n"
for l in readme:
    if '```' in l and not in_code:
        in_code = True
        continue
    if '```' not in l and in_code:
        code_lines.append(l)
        continue
    elif '```' in l and in_code:
        in_code = False
        # don't need to add newlines
        _index_html_template_ += f"<pre><code>{newline.join(code_lines)}</code></pre>\n"
        code_lines = []
        continue

    if l.lstrip().startswith('-'):
        in_list = True
        list_lines.append(l)
        continue
    if in_list and not l.lstrip().startswith('-'):
        in_list = False
        # don't need to add newlines
        _index_html_template_ += f"<pre><code>\n{newline.join(list_lines)}</code></pre>\n"
        list_lines = []


    # headings
    if '\_' in l:
        l.replace('\_', '_')
    # insert the actual text as a paragraph
    if l.startswith("#"):
        poundcount = len(l.split()[0])
        if poundcount <= 6:
            _index_html_template_ += f"<h{poundcount}>{l.replace('#'*poundcount, '')}</h{poundcount}>\n"
        else:
            print("warning poundcount > 6!, defaulting to paragraph")
            _index_html_template_ += "<p>" + l + "</p>\n"
    else:
        _index_html_template_ += "<p>" + l + "</p>\n"

# html closing tag
_index_html_template_ += "</body></html>\n"

print(_index_html_template_)

with open(f'{working_dir}/index.html', 'w') as f:
    f.write(_index_html_template_)

# write the overall results
with open(_multi_step_result_, 'w') as f:
    if len(failures) == 0 or debug:
        overall_result = _S_passed_
    else:
        overall_result = _S_failed_
    r = {
        _StateKey_  : overall_result,
        _ResultsKey_ : {},
        _LogKey_ : 'Done!',
    }

    json.dump(r, f, sort_keys=True, indent=2)


benchmark.save_variables()  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
