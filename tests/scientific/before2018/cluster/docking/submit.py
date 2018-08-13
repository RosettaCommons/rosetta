#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import os, sys, commands


class NT:  # named tuple
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = '|'
        for i in dir(self):
            if not i.startswith('__'): r += '%s --> %s, ' % (i, getattr(self, i))
        return r[:-2]+'|'


print 'Running submit.py script for cluster loop scientific test...'

m_vars = eval( file('_arguments.py').read() )

targets = {}
for line in file('input/sources.txt'):
    s = line.split()
    if s: targets[ s[0] ] = NT(partners=s[1], pH=s[2])

# Prepack...
prepack_template = '%(bin)s/docking_prepack_protocol.%(binext)s -database %(database)s \
-s input/%(target)s.u.opt.pdb -out::path:pdb output/%(target)s/prepack/ \
-partners %(partners)s -dock_ppk -ex1 -ex2aro -unboundrot input/%(target)s.u.opt.pdb -out:level 200'

for t in targets:
    os.mkdir('output/'+t);  os.mkdir('output/'+t+'/prepack')
    m_vars['target'] = t
    m_vars['partners'] = targets[t].partners

    commandline = prepack_template % m_vars

    res, output = commands.getstatusoutput(commandline)
    print output
    if res: sys.exit(res)


T_condor_job_header = '''
universe     = vanilla
Notify_user  =
notification = Error
Log          = .condorscript.log
Executable   = %(bin)s/docking_protocol.%(binext)s

Requirements = ( Memory > 256)
GetEnv       = True

##on_exit_remove = (ExitBySignal == False)

'''


T_condor_job_target = '''
## now looprelax jobs
Error   = output/%(target)s.logerr
Output  = output/%(target)s.logout
arguments = -database %(database)s -in::file::s output/%(target)s/prepack/%(target)s.u.opt_0001.pdb -nstruct 1000 \
            -dock_pert 3 8 -spin -out::path:pdb output/%(target)s/ -partners %(partners)s -ex1 -ex2aro \
            -unboundrot input/%(target)s.u.opt.pdb -native input/%(target)s.b.opt.pdb \
            -run:multiple_processes_writing_to_one_directory -out:file:scorefile ./output/%(target)s/%(target)s.fasc \
            -random_delay 5 -seed_offset $(Process)

priority = -10

queue 1

'''
#-nstruct 500

condor = T_condor_job_header % m_vars

for t in targets:
    #os.mkdir('output/'+t)
    m_vars['target'] = t
    m_vars['partners'] = targets[t].partners
    condor += T_condor_job_target % m_vars

f = file('condor', 'w');  f.write(condor);  f.close()
