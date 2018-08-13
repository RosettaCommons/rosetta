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
Executable   = %(bin)s/rosetta_scripts.%(binext)s

Requirements = ( Memory > 512)
GetEnv       = True

'''

T_condor_job_target = '''
Error   = output/%(target)s-r.logerr
Output  = output/%(target)s-r.logout

arguments = @flags.txt -nstruct 100 -database %(database)s \\
            -in:file:s input/%(target)s.pdb.gz -in:file:native native/%(target)s.pdb.gz \\
            -packing:unboundrot unbound_from_kwk/%(target)s.pdb.gz \\
						-out:path:pdb output -out:file:atom_tree_diff %(target)s-r-atom_tree_diff.out \\
						-parser:protocol rotate.xml

priority = -10
queue 1
'''


T_condor_job2_target = '''
Error   = output/%(target)s-rt-%(n)s.logerr
Output  = output/%(target)s-rt-%(n)s.logout

arguments = @flags.txt -nstruct 100 -database %(database)s \\
            -in:file:s input/%(target)s.pdb.gz -in:file:native native/%(target)s.pdb.gz \\
            -packing:unboundrot unbound_from_kwk/%(target)s.pdb.gz -run:seed_offset %(n)s \\
						-out:path:pdb output -out:file:atom_tree_diff %(target)s-rt-atom_tree_diff.%(n)s \\
						-parser:protocol translate_rotate.xml

priority = -10
queue 1
'''

print 'Running submit.py script for cluster ligand_dock_script scientific test...'

m_vars = eval( file('_arguments.py').read() )

condor = T_condor_job_header % m_vars

os.mkdir('output')

for t in targets:
    #os.mkdir('output/'+t)
    m_vars['target'] = t
    condor += T_condor_job_target % m_vars
    for n in range(6):
        m_vars['n'] = n
        condor += T_condor_job2_target % m_vars

f = file('condor', 'w')
f.write(condor)
f.close()
