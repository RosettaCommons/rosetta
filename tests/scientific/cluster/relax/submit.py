#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# --> relax <-- scientific test


import os

from targets import targets

T_condor_job_header = '''
universe     = vanilla
Notify_user  =
notification = Error
Log          = .condorscript.log
Executable   = %(bin)s/relax.%(binext)s

Requirements = ( Memory > 256)
GetEnv       = True

'''


T_condor_job_target = '''
## now looprelax jobs
Error   = output/%(target)s/logerr
Output  = output/%(target)s/logout
arguments = -database %(database)s -run:protocol relax -relax:cycle_ratio 1.0 \
            -l input/centroid_abinitio/%(target)s/list -out:file:silent output/%(target)s/%(target)s.out \
            -in:file:fullatom -native input/fragments/%(target)s/boinc_rb1_%(target)s.pdb \
            -random_delay 5 -silent_decoytime -mute all -seed_offset $(Process)

priority = -10

queue 1

'''


print 'Running submit.py script for cluster loop scientific test...'


m_vars = eval( file('_arguments.py').read() )

condor = T_condor_job_header % m_vars

for t in targets:
    os.mkdir('output/'+t)
    m_vars['target'] = t
    condor += T_condor_job_target % m_vars

f = file('condor', 'w');  f.write(condor);  f.close()
