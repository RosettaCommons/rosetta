#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import os

from targets import targets


T_condor_job_header = '''
universe     = vanilla
Notify_user  =
notification = Error
Log          = .condorscript.log
Executable   = %(bin)s/loopmodel.%(binext)s

Requirements = ( Memory > 256)
GetEnv       = True

##on_exit_remove = (ExitBySignal == False)

'''


T_condor_job_target = '''
## now looprelax jobs
Error   = output/%(target)s.logerr
Output  = output/%(target)s.logout
arguments = -database %(database)s -loops::frag_sizes 9 3 1 -out:nstruct 25 -score:weights score13_env_hb \
            -loops::frag_files input/fragments/%(target)s/frag9.gz ./input/fragments/%(target)s/frag3.gz none \
            -in::file::s input/fragments/%(target)s/%(target)s.pdb \
            -loops::loop_file input/loopdefs/%(target)s.loopfile \
            -in:file:native input/fragments/%(target)s/%(target)s.pdb \
            -psipred_ss2 input/fragments/%(target)s/psipred_ss2 -loops::extended \
            -out:file:silent output/%(target)s.out -out::prefix %(target)s \
            -loops:remodel quick_ccd -loops:refine no -out:file:silent_decoytime -mute all \
            -random_delay 5 -seed_offset $(Process)

priority = -10

queue 20

'''
#-nstruct 500


print 'Running submit.py script for cluster loop scientific test...'


m_vars = eval( file('_arguments.py').read() )

condor = T_condor_job_header % m_vars

for t in targets:
    #os.mkdir('output/'+t)
    m_vars['target'] = t
    condor += T_condor_job_target % m_vars

f = file('condor', 'w');  f.write(condor);  f.close()
