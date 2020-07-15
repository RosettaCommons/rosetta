#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  peptide_pnear_vs_ic50/9.finalize.py
## @brief this script is part of peptide_pnear_vs_ic50 scientific test
## @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).


import json, subprocess
import benchmark

from benchmark import *

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

# Read readme
readme = subprocess.getoutput( "cat readme.md" ).splitlines()

# Read previous analysis
all_peptide_analyses = []
all_failures = []
suffices = [ "1A", "1B", "1C", "1D", "1E", "1F", "1G" ]
overall_pass = True
at_least_one_peptide_fails = False
for suffix in  suffices:
    curanalysis = subprocess.getoutput( "cat result_" + suffix + ".txt" ).splitlines()
    curfailures = []
    for line in curanalysis:
        linesplit = line.split()
        if  len(linesplit) == 0 :
            continue
        if linesplit[ len(linesplit) - 1 ] == "NO" and linesplit[0] != "OVERALL" :
            overall_pass = False
            at_least_one_peptide_fails = True
            curfailures.append(line)
    all_peptide_analyses.append(curanalysis)
    all_failures.append(curfailures)

# Get R^2 value from fit and confirm that it's greater than 0.85.
Rsq_from_fit = bundle[2]
Rsq_passes = True
rsq_string = "R-squared value from fitting greater than 0.85?\tYES"
if Rsq_from_fit < 0.85 :
    Rsq_passes = False
    overall_pass = False
    rsq_string = "R-squared value from fitting greater than 0.85?\tNO"

# Build up html from readme, starting with failures:
_index_html_template_ = "<html>\n"
_index_html_template_ += "<H2>Scientific test: " + testname + "</H2>\n"
_index_html_template_ += "<h3>FAILURES</h3>\n<p>\n"
if overall_pass == False:
    if at_least_one_peptide_fails == True :
        _index_html_template_ += "The following peptides showed failures:" + "<br/>\n"
        for i in range(len(all_failures)) :
            curfailures = all_failures[i]
            if len(curfailures) == 0 :
                continue
            _index_html_template_ += "&nbsp;&nbsp;&nbsp;&nbsp;Peptide NDM1i_" + suffices[i] + ":<br/>\n"
            for failure in curfailures :
                _index_html_template_ += "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" + failure + "<br/>\n"
    if Rsq_passes == False :
        _index_html_template_ += rsq_string
else:
    _index_html_template_ += "None." + "<br/>\n"

# add pass/fail summary to html
_index_html_template_ += "<h3>RESULTS SUMMARY</h3>\n<p>\n"

# add image for results
_index_html_template_ += 'Analysis of correlation of predicted folding propensity with experimentally-measured inhibition values (which should be linear):<br/>\n<img src="fitcurve.png" width="600" alt="Correlation of delta-G of folding (computed from the funnels above) with IC50 values (measured experimentally).  This should be roughly linear." style="max-width: 100%" /><br/><br/>\n'
_index_html_template_ += 'Folding funnels:<br/>\n<img src="e_vs_rmsd_plots.png" width="800" alt="Energy landscapes computed for peptide to designed structure (left) or to lowest-energy structure sampled (right).  These should look funnel-shaped." style="max-width: 100%" /><br/><br/>\n'

for i in range(len(suffices)) :
    suffix = suffices[i]
    analysis = all_peptide_analyses[i]
    _index_html_template_ += "Peptide NDM1i_" + suffix + ":<br/>\n"
    firstline = True
    for line in analysis:
        # Skip the first line, which redundantly lists the peptide:
        if firstline == True:
            firstline = False
            continue
        _index_html_template_ += "&nbsp;&nbsp;&nbsp;&nbsp;" + line + "<br/>\n"
    _index_html_template_ += "<br/>\n"
_index_html_template_ += "Correlation analysis:<br/>\n&nbsp;&nbsp;&nbsp;&nbsp;" + rsq_string + "<br/><br/>\n"

# add text from readme
for l in readme:

    # headings
    if l.startswith( "## " ):
        _index_html_template_ += "<h3>" + l.replace( "## ", "" ) + "</h3>\n"
    
    # ignore the description
    elif l.startswith( "#### " ):
        _index_html_template_ += "<p><i>" + l.replace( "#### ", "" ) + "</i></p>\n"

    # correct the image

    elif l == "![Improvements from talaris2013 through ref2015](inputs/Mulligan_2020_SuppFig_comparing_sfxns.png)" :
        _index_html_template_ += '<img src="inputs/Mulligan_2020_SuppFig_comparing_sfxns.png" width="800" alt="Improvements from talaris2013 through ref2015." style="max-width: 100%" /><br/><br/>\n'
        
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
