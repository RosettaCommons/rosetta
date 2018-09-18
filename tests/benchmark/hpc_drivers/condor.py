# -*- coding: utf-8 -*-
# :noTabs=true:

import os, sys, time, stat


import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) + '/base.py')  # A bit of Python magic here, what we trying to say is this: from base import *, but path to base is calculated from our source location  # from base import HPC_Driver, execute, NT


T_condor_job_template = '''
universe     = vanilla
Notify_user  =
notification = Error
Log          = {log_dir}/.hpc.{name}.condor.log
Executable   = {executable}

periodic_remove = JobStatus == 5

request_memory = {memory}

GetEnv       = True

# Target: {name}
Output  = {log_dir}/.hpc.{name}.output.$(Process).log
Error   = {log_dir}/.hpc.{name}.errors.$(Process).log

arguments = {arguments}

InitialDir = {working_dir}

# Removing jobs that run for more then specified number of seconds
periodic_remove = (RemoteWallClockTime - CumulativeSuspensionTime) > {run_time}

queue {jobs_to_queue}
'''
#Requirements =  (Memory > 256)  &&  (Arch == "X86_64")
#Requirements = Arch == "X86_64"
#Requirements = (Memory > %(memory)s) %(requirements)s

# requirelements for a new nodes:
#Requirements = (Machine > "rosetta012.graylab.jhu.edu")

# requirelements for old nodes:
# Requirements = (Machine <= "rosetta012.graylab.jhu.edu")


T_condor_target_template = '''
'''
#priority = -10
#T_condor_job_target = '''


class Condor_HPC_Driver(HPC_Driver):
    def get_condor_accumulated_usage(self, user='$USER@'):
        # Expected output: sergey@UT64 0.50 0.50 1.00 0 196.86 12/26/2010 23:30  9/27/2011 23:55
        for _ in range(8):
            try:
                o = execute('', 'condor_userprio -all -allusers | grep {}'.format(user), return_='output', terminate_on_failure=False).split()
                if len(o) >= 6: return max(1.0, float( o[5] ) )
                else: return 0.0
            except ValueError as _: pass

            time.sleep(32)

        return 0.0



    # def complete(self, condor_job_id):
    #     ''' Return job completion status. Note that single hpc_job may contatin inner list of individual HPC jobs, True should be return if they all run in to completion.
    #     '''

    #     execute('Releasing condor jobs...', 'condor_release $USER', return_='tuple')

    #     s = execute('', 'condor_q $USER | grep $USER | grep {}'.format(condor_job_id), return_='output', terminate_on_failure=False).replace(' ', '').replace('\n', '')
    #     if s: return False

    #         # #setDaemonStatusAndPing('[Job #%s] Running... %s condor job(s) in queue...' % (self.id, len(s.split('\n') ) ) )
    #         # n_jobs = len(s.split('\n'))
    #         # s, o = execute('', 'condor_userprio -all | grep $USER@', return_='tuple')
    #         # if s == 0:
    #         #     jobs_running = o.split()
    #         #     jobs_running = 'XX' if len(jobs_running) < 4 else jobs_running[4]
    #         #     self.set_daemon_message("Waiting for condor to finish HPC jobs... [{} jobs in HPC-Queue, {} CPU's used]".format(n_jobs, jobs_running) )
    #         #     print "{} condor jobs in queue... Sleeping 32s...    \r".format(n_jobs),
    #         # sys.stdout.flush()
    #         # time.sleep(32)
    #     else:

    #         #self.tracer('Waiting for condor to finish the jobs... DONE')
    #         self.jobs.remove(condor_job_id)
    #         self.cpu_usage += self.get_condor_accumulated_usage()
    #         return True  # jobs already finished, we return empty list to prevent double counting of cpu_usage


    def complete(self, condor_job_id):
        ''' Return True if job with given id is complete
        '''
        execute('Releasing condor jobs...', 'condor_release $USER', return_='tuple', silent=True)

        s = execute('', f'condor_q $USER | grep $USER | grep {condor_job_id}', return_='output', terminate_on_failure=False, silent=True)
        if s: return False
        else:
            #self.tracer('Waiting for condor to finish the jobs... DONE')
            self.jobs.remove(condor_job_id)
            self.cpu_usage += self.get_condor_accumulated_usage()
            return True  # jobs already finished, we return empty list to prevent double counting of cpu_usage


    def cancel_job(self, condor_job_id):
        execute(f'Condor_HPC_Driver.canceling job {condor_job_id}...', f'condor_rm {condor_job_id}', terminate_on_failure=False)




    def submit_hpc_job(self, name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory=512, time=12, block=True):
        self.cpu_usage -= self.get_condor_accumulated_usage()

        # creating shell wrapper in order to reliably capture output
        #execute_sh = os.path.abspath( self.working_dir + '/.hpc.execute.{}.sh'.format(name) )
        #with open(execute_sh, 'w') as f: f.write('#!/bin/bash\n$*\n');  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)

        arguments = arguments.format(process='$(Process)')
        run_time=int(time*60*60)

        #jobs = []
        #???? p = dict(j, , arguments=j['arguments'].format(process='$(Process)') )

        #execute_sh = os.path.abspath(self.working_dir + '/.hpc.execute.{}.sh'.format(j['target']))
        #with open(execute_sh, 'w') as f: f.write('#!/bin/bash\n{executable} {arguments}\n'.format(**j));  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)

        condor_file = working_dir + '/.hpc.{}.condor'.format(name)
        condor_spec = T_condor_job_template.format(name=name, executable=executable, arguments=arguments, working_dir=working_dir, jobs_to_queue=jobs_to_queue, log_dir=log_dir,
                                                   memory=memory, process='$(Process)', run_time=run_time)
                                                   #requirements=self.config.get('condor', 'requirements'),
                                                   #execute_sh=execute_sh)

        with open(condor_file, 'w') as f: f.write(condor_spec)
        condor_job_id = int( execute('Submitting jobs to condor...', 'cd {} && condor_submit {}'.format(self.working_dir, condor_file),
                                     tracer=self.tracer, return_='output').split()[-1][:-1] )

          #jobs.append( NT(condor_file=condor_file, condor_job_id=condor_job_id, ) )

        self.jobs.append(condor_job_id)

        if block:
            self.wait_until_complete( [condor_job_id] )
            return None

        else: return condor_job_id
