#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

import os

from targets import targets

T_condor_job_header = '''
universe     = vanilla
Notify_user  =
notification = Error
Log          = .condorscript.log
Executable   = %(bin)s/docking_protocol.%(binext)s

Requirements = ( Memory > 512)
GetEnv       = True

'''

T_condor_job_target = '''
Error   = output/%(target)s.logerr
Output  = output/%(target)s.logout

arguments = -database %(database)s \
            -s input/%(target)s.pdb -out:file:scorefile %(target)s -out:path:score output -out:path:pdb output/decoys \
            -multiple_processes_writing_to_one_directory  \
            -dock_pert 3 8 -spin -ex1 -ex2aro -nstruct 1000 -out:prefix sci_bench -out:file:fullatom -mute core

priority = -10
queue 50

'''


print 'Running submit.py script for cluster docking-2 scientific test...'

m_vars = eval( file('_arguments.py').read() )

condor = T_condor_job_header % m_vars

os.mkdir('output')

for t in targets:
    m_vars['target'] = t
    condor += T_condor_job_target % m_vars

f = file('condor', 'w');  f.write(condor);  f.close()
