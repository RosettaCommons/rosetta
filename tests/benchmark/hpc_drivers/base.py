# -*- coding: utf-8 -*-
# :noTabs=true:

import os, commands, stat
import time as time_module


class NT:  # named tuple
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = 'NT: |'
        for i in dir(self):
            if not i.startswith('__') and not isinstance(getattr(self, i), types.MethodType): r += '{} --> {}, '.format(i, getattr(self, i))
        return r[:-2]+'|'



class HPC_Exception(Exception):
    def __init__(self, value): self.value = value
    def __str__(self): return self.value



def execute(message, commandline, return_=False, until_successes=False, tracer=lambda x:None, terminate_on_failure=True):
    tracer(message);  tracer(commandline)
    while True:
        (res, output) = commands.getstatusoutput(commandline)
        tracer(output)

        if res and until_successes: pass  # Thats right - redability COUNT!
        else: break

        print "Error while executing %s: %s\n" % (message, output)
        print "Sleeping 60s... then I will retry..."
        time.sleep(60)

    if return_ == 'tuple': return(res, output)

    if res and terminate_on_failure:
        tracer("\nEncounter error while executing: " + commandline)
        if return_==True: return True
        else: raise HPC_Exception("\nEncounter error while executing: " + commandline + '\n' + output)

    if return_ == 'output': return output
    else: return False


def Sleep(time_, message, dict_={}):
    ''' Fancy sleep function '''
    len_ = 0
    for i in range(time_, 0, -1):
        #print "Waiting for a new revision:%s... Sleeping...%d     \r" % (sc.revision, i),
        msg = message.format( **dict(dict_, time_left=i) )
        print msg,
        len_ = max(len_, len(msg))
        sys.stdout.flush()
        time_module.sleep(1)

    print ' '*len_ + '\r',  # erazing sleep message


# Abstract class for HPC job submission
class HPC_Driver:
    def __init__(self, working_dir, config, tracer=lambda x:None, set_daemon_message=lambda x:None):
        self.working_dir = working_dir
        self.config = config
        self.cpu_usage  = 0.0  # cummulative cpu usage in hours
        self.tracer     = tracer
        self.set_daemon_message = set_daemon_message

        self.cpu_count = self.config.getint('DEFAULT', 'cpu_count')

        self.jobs = []  # list of all jobs currently running by this driver, Job class is driver depended, could be just int or something more complex


    def execute(self, executable, arguments, working_dir, log_dir=None, name='_no_name_', memory=256, time=24, shell_wrapper=False, block=True):
        ''' Execute given command line on HPC cluster, must accumulate cpu hours in self.cpu_usage '''
        if log_dir==None: log_dir=self.working_dir

        if shell_wrapper:
            shell_wrapper_sh = os.path.abspath(self.working_dir + '/hpc.{}.shell_wrapper.sh'.format(name))
            with file(shell_wrapper_sh, 'w') as f: f.write('#!/bin/bash\n{} {}\n'.format(executable, arguments));  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)
            executable, arguments = shell_wrapper_sh, ''

        return self.submit_hpc_job(name=name, executable=executable, arguments=arguments, working_dir=working_dir, log_dir=log_dir, jobs_to_queue=1, memory=memory, time=time, block=block)


    def submit_hpc_job(self, name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory=512, time=12, block=True):
        must_be_implemented_in_inherited_classes


    def complete(self, job):
        ''' Return job completion status. Note that single hpc_job may contatin inner list of individual HPC jobs, True should be return if they all run in to completion.
        '''
        must_be_implemented_in_inherited_classes


    def wait_until_complete(self, jobs=None):
        ''' Helper function, wait until given jobs list is finished, if no argument is given waits until all jobs known by driver is finished '''
        jobs = jobs if jobs else self.jobs

        while jobs:
            for j in jobs:
                if self.complete(j): jobs.remove(j)

            if jobs:
                #total_cpu_queued  = sum( [j.jobs_queued  for j in jobs] )
                #total_cpu_running = sum( [j.cpu_running for j in jobs] )
                #self.set_daemon_message("Waiting for HPC job(s) to finish... [{} process(es) in queue, {} process(es) running]".format(total_cpu_queued, total_cpu_running) )
                #self.tracer("Waiting for HPC job(s) [{} process(es) in queue, {} process(es) running]...  \r".format(total_cpu_queued, total_cpu_running), end='')
                #print "Waiting for {} HPC jobs to finish... [{} jobs in queue, {} jobs running]... Sleeping 32s...     \r".format(total_cpu_queued, cpu_queued+cpu_running, cpu_running),

                self.set_daemon_message("Waiting for HPC {} job(s) to finish...".format( len(jobs) ) )
                #self.tracer("Waiting for HPC {} job(s) to finish...".format( len(jobs) ) )

                sys.stdout.flush()
                Sleep(64, '"Waiting for HPC {n_jobs} job(s) to finish, sleeping {time_left}s    \r', dict(n_jobs=len(jobs)))
                #time_module.sleep(64*1)
