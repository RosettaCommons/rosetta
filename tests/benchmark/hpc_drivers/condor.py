# -*- coding: utf-8 -*-
# :noTabs=true:

import os, sys, time, stat


import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) + '/base.py')  # A bit of Python magic here, what we trying to say is this: from base import *, but path to base is calculated from our source location  # from base import HPC_Driver, execute, NT


T_condor_job_template = '''
universe     = vanilla
Notify_user  =
notification = Error
Log          = {log_dir}/.hpc.{target}.condor.log
Executable   = {execute_sh}

periodic_remove = JobStatus == 5

request_memory = {memory}

GetEnv       = True


# Target: {target}
Output  = {log_dir}/.hpc.{target}.output.$(Process).log
Error   = {log_dir}/.hpc.{target}.errors.$(Process).log

arguments = {executable} {arguments}

InitialDir = {initial_dir}

# Removing jobs that run for more then specified number of seconds
periodic_remove = (RemoteWallClockTime - CumulativeSuspensionTime) > {run_time}

queue {queue}
'''
#Requirements =  (Memory > 256)  &&  (Arch == "X86_64")
#Requirements = Arch == "X86_64"
#Requirements = (Memory > %(memory)s) %(requirements)s


T_condor_target_template = '''
'''
#priority = -10
#T_condor_job_target = '''


class Condor_HPC_Driver(HPC_Driver):
    def get_condor_accumulated_usage(self, user='$USER@'):
        # Expected output: sergey@UT64 0.50 0.50 1.00 0 196.86 12/26/2010 23:30  9/27/2011 23:55
        o = execute('', 'condor_userprio -all -allusers | grep {}'.format(user), return_='output', terminate_on_failure=False).split()
        if len(o) >= 6: return max(1.0, float( o[5] ) )
        else: return 0.0


    def execute_hpc_jobs(self, hpc_jobs):
        self.cpu_usage -= self.get_condor_accumulated_usage()

        # creating shell wrapper in order to reliably capture output
        execute_sh = os.path.abspath(self.working_dir + '/.hpc.execute.sh')
        with file(execute_sh, 'w') as f: f.write('#!/bin/bash\n$*\n');  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)


        jobs = []
        for j in hpc_jobs:
            p = dict(j, run_time=int(j['run_time']*60*60), arguments=j['arguments'].format(process='$(Process)') )

            #execute_sh = os.path.abspath(self.working_dir + '/.hpc.execute.{}.sh'.format(j['target']))
            #with file(execute_sh, 'w') as f: f.write('#!/bin/bash\n{executable} {arguments}\n'.format(**j));  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)

            condor_file = self.working_dir + '/.hpc.{}.condor'.format(p['target'])
            condor_spec = T_condor_job_template.format(process='$(Process)', log_dir=self.working_dir,
                                                       #requirements=self.config.get('condor', 'requirements'),
                                                       execute_sh=execute_sh, **p)


            with file(condor_file, 'w') as f: f.write(condor_spec)
            condor_job_id = int( execute('Submitting jobs to condor...', 'cd %s && condor_submit %s' % (self.working_dir, condor_file),
                                         tracer=self.tracer, return_='output').split()[-1][:-1] )

            jobs.append(condor_job_id)  #jobs.append( NT(condor_file=condor_file, condor_job_id=condor_job_id, ) )

        if True:  # if block:
            self.tracer('Waiting for condor to finish the jobs...')
            while True:
                execute('Releasing condor jobs...', 'condor_release $USER', return_='tuple')
                s = ''
                for j_id in jobs: s += execute('', 'condor_q $USER | grep $USER | grep {}'.format(j_id), return_='output', terminate_on_failure=False)
                if s:
                    #setDaemonStatusAndPing('[Job #%s] Running... %s condor job(s) in queue...' % (self.id, len(s.split('\n') ) ) )
                    n_jobs = len(s.split('\n'))
                    s, o = execute('', 'condor_userprio -all | grep $USER@', return_='tuple')
                    if s == 0:
                        jobs_running = o.split()
                        jobs_running = 'XX' if len(jobs_running) < 4 else jobs_running[4]
                        self.set_daemon_message("Waiting for condor to finish HPC jobs... [{} jobs in HPC-Queue, {} CPU's used]".format(n_jobs, jobs_running) )
                        print "{} condor jobs in queue... Sleeping 32s...    \r".format(n_jobs),
                    sys.stdout.flush()
                    time.sleep(32)
                else: break
            self.tracer('Waiting for condor to finish the jobs... DONE')
            self.cpu_usage += self.get_condor_accumulated_usage()
            return []  # jobs already finished, we return empty list to prevent double counting of cpu_usage

        #self.jobs.append(condor_job_ids)
        #return condor_job_ids
