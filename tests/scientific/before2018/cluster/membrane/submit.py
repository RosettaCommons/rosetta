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
Executable   = %(bin)s/membrane_abinitio2.%(binext)s

Requirements = ( Memory > 256)
GetEnv       = True

'''

T_condor_job_target = '''
Error   = output/%(target)s.logerr
Output  = output/%(target)s.logout

arguments = -database %(database)s -abinitio:membrane -mute core.scoring.MembranePotential \
            -nstruct 5 \
            -frag9 input/%(target)s/aa%(target)s_09_05.200_v1_3 \
            -frag3 input/%(target)s/aa%(target)s_03_05.200_v1_3 \
            -in::file::native input/%(target)s/%(target)s.pdb \
            -in::file::spanfile input/%(target)s/%(target)s.span \
            -in::file::lipofile input/%(target)s/%(target)s.lips4 \
            -membrane:no_interpolate_Mpair \
            -membrane:Menv_penalties \
            -score:find_neighbors_3dgrid \
            -random_delay 7 -seed_offset $(Process) \
            -multiple_processes_writing_to_one_directory \
            -out:file:silent output/%(target)s.out

priority = -10

queue 80

'''



print 'Running submit.py script for membrane abinitio scientific test...'

m_vars = eval( file('_arguments.py').read() )

condor = T_condor_job_header % m_vars

os.mkdir('output')

for t in targets:
    #os.mkdir('output/'+t)
    m_vars['target'] = t
    #m_vars['target_short'] = t[:-1]
    condor += T_condor_job_target % m_vars

f = file('condor', 'w');  f.write(condor);  f.close()
