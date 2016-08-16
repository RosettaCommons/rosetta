#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import os.path, sys, commands

#import yaml
import json

class CD:
    def __init__(self, **entries): self.__dict__.update(entries)


results = {}

res, output = commands.getstatusoutput('./post_process.sh')
print output
if res:
    print 'Abnormal termination in: ./post_process.sh, code:', res
    sys.exit(res)

lines = file('./output/results.txt').read().split('\n')

mean_n5_index = lines[0].split(';').index('mean_n5')
sd_n5_index   = lines[0].split(';').index('sd_n5')

num_of_passed_targets = 0
for l in lines[1:]:
    fields = l.split(';')
    if len(fields) > max(mean_n5_index, sd_n5_index):
        if float(fields[mean_n5_index]) > 3.0 : num_of_passed_targets += 1
        if fields[0].split('/') > 2:
            target = fields[0].split('/')[2]
            results[ target + '.mean_n5' ] = float(fields[mean_n5_index])
            results[ target + '.sd_n5' ]   = float(fields[sd_n5_index])

results['num_of_passed_targets'] = num_of_passed_targets
results['_isTestPassed'] = num_of_passed_targets > 40


f = file('.results.yaml', 'w');  f.write( json.dumps(results) );  f.close()
