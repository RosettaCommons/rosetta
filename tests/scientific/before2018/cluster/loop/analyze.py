#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import os.path, commands

import yaml

from targets import targets

class CD:
    def __init__(self, **entries): self.__dict__.update(entries)


results = {}

for t in targets:
    fname = 'output/%s.out' % t
    if os.path.isfile(fname):
        (res, output) = commands.getstatusoutput('cd output && cat %s.out | grep SCORE:' % t)
        lines = output.split('\n')
        legend = lines[0].split()
        keys = ['score', 'rms']
        for k in keys:
            if k not in legend: continue  # we check if all needed info are here. If not - just skip it...

        rms_key = legend.index('rms')
        score_key = legend.index('score')

        lowest_rmsd = 1e100
        rmsd_of_the_lowest_energy_structure = 1e100
        lowest_score = 1e100

        for l in lines[1:]:
            l_split = l.split()
            #if (rms_key not in l_split) or (score_key not in l_split): continue

            rms = float( l_split[rms_key] )
            score = float( l_split[score_key] )

            if lowest_rmsd > rms : lowest_rmsd = rms

            if lowest_score > score :
                lowest_score = score
                rmsd_of_the_lowest_energy_structure = rms

        results[t] = dict(lowest_rmsd=lowest_rmsd, lowest_score=lowest_score, rmsd_of_the_lowest_energy_structure=rmsd_of_the_lowest_energy_structure)

# Now, compare results with pre-defined thresholds
thresholds = yaml.load( file('ScientificCluster.Loop.thresholds').read() )
num_of_passed_targets = 0
for k in thresholds:
    for sk in thresholds[k]:
        if k in results  and  sk in results[k]  and  results[k][sk] < thresholds[k][sk]: num_of_passed_targets += 1

results['num_of_passed_targets'] = num_of_passed_targets
results['_isTestPassed'] = num_of_passed_targets > 34


f = file('.results.yaml', 'w');  f.write( yaml.dump(results) );  f.close()
