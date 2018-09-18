# -*- coding: utf-8 -*-
# :noTabs=true:

import time as time_module
import codecs
import signal

import sys, imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) + '/base.py')  # A bit of Python magic here, what we trying to say is this: from base import *, but path to base is calculated from our source location  # from base import HPC_Driver, execute, NT




class MultiCore_HPC_Driver(HPC_Driver):

    class JobID:
        def __init__(self, pids=None):
            self.pids = pids if pids else []


        def __bool__(self): return bool(self.pids)


        def __len__(self): return len(self.pids)


        def add_pid(self, pid): self.pids.append(pid)


        def remove_completed_pids(self):
            for pid in self.pids[:]:
                try:
                    r = os.waitpid(pid, os.WNOHANG)
                    if r == (pid, 0): self.pids.remove(pid)  # process have ended without error
                    elif r[0] == pid :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                        self.cancel_job()
                        print('Some of the HPC jobs terminated abnormally!')
                        sys.exit(1)

                except ChildProcessError: self.pids.remove(pid)


        def cancel(self):
            for pid in self.pids:
                try:
                    os.killpg(os.getpgid(pid), signal.SIGKILL)
                except ChildProcessError: pass

            self.pids = []



    def __init__(self, *args, **kwds):
        HPC_Driver.__init__(self, *args, **kwds)
        #print(f'MultiCore_HPC_Driver: cpu_count: {self.cpu_count}')


    def remove_completed_jobs(self):
        for job in self.jobs[:]:  # Need to make a copy so we don't modify a list we're iterating over
            job.remove_completed_pids()
            if not job: self.jobs.remove(job)


    @property
    def process_count(self):
        ''' return number of processes that currently ran by this driver instance
        '''
        return sum( map(len, self.jobs) )


    def submit_hpc_job(self, name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory=512, time=12, block=True):
        cpu_usage = -time_module.time()/60./60.

        def mfork():
            ''' Check if number of child process is below cpu_count. And if it is - fork the new pocees and return its pid.
            '''
            while self.process_count >= self.cpu_count:
                self.remove_completed_jobs()
                if self.process_count >= self.cpu_count: time_module.sleep(.5)

            sys.stdout.flush()
            pid = os.fork()
            # appending at caller level insted  if pid: self.jobs.append(pid) # We are parent!
            return pid

        current_job = self.JobID()
        process = 0
        for i in range(jobs_to_queue):

            pid = mfork()
            if not pid: # we are child process
                command_line = 'cd {} && {} {}'.format(working_dir, executable, arguments.format(process=process) )
                log = execute('Running job {}.{}...'.format(name, i), command_line, tracer=self.tracer, return_='output')
                with codecs.open(log_dir+'/hpc.{name}.{i}.log'.format(**vars()), 'w', encoding='utf-8', errors='replace') as f: f.write(command_line+'\n'+log)
                sys.exit(0)
            else: # we are parent!
                current_job.add_pid(pid)
                # Need to potentially re-add to list, as remove_completed_jobs() might trim it.
                if current_job not in self.jobs: self.jobs.append(current_job)

            process += 1

        if block:
            #for p in all_queued_jobs: os.waitpid(p, 0)  # waiting for all child process to termintate...

            self.wait_until_complete(current_job)
            self.remove_completed_jobs()

            cpu_usage += time_module.time()/60./60.
            self.cpu_usage += cpu_usage * jobs_to_queue  # approximation...

            current_job = self.JobID()

        return current_job


    def complete(self, job_id):
        ''' Return job completion status. Return True if job completed and False otherwise
        '''
        self.remove_completed_jobs()
        return job_id not in self.jobs


    def cancel_job(self, job):
        job.cancel();
        if job in self.jobs:
            self.jobs.remove(job)


    def __repr__(self):
        return 'MultiCore_HPC_Driver<>'
