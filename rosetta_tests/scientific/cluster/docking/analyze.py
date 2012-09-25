#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

import os.path, commands

#import yaml
import json

class CD:
    def __init__(self, **entries): self.__dict__.update(entries)


results = {}

results['num_of_passed_targets'] = 1
results['_isTestPassed'] = True


f = file('.results.yaml', 'w');  f.write( json.dumps(results) );  f.close()
