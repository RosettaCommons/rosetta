# -*- coding: utf-8 -*-
# :noTabs=true:

import os, commands, stat

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



def execute(message, commandline, return_=False, untilSuccesses=False, tracer=lambda x:None, terminate_on_failure=True):
    tracer(message);  tracer(commandline)
    while True:
        (res, output) = commands.getstatusoutput(commandline)
        tracer(output)

        if res and untilSuccesses: pass  # Thats right - redability COUNT!
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


# Abstract class for HPC job submission
class HPC_Driver:
    def __init__(self, working_dir, config, tracer=lambda x:None, set_daemon_message=lambda x:None):
        self.working_dir = working_dir
        self.config = config
        self.cpu_usage  = 0.0  # cummulative cpu usage in hours
        self.tracer     = tracer
        self.set_daemon_message = set_daemon_message

        self.jobs = []  # list of all jobs currently running by this driver, Job class is driver depended, could be just int or something more complex


    def execute(self, command_line, initial_dir, target, memory=256, run_time=24, shell_wrapper=True):
        ''' Execute given command line on HPC cluster, must accumulate cpu hours in self.cpu_usage '''

        if shell_wrapper:
            shell_wrapper_sh = os.path.abspath(self.working_dir + '/hpc.{}.shell_wrapper.sh'.format(target))
            with file(shell_wrapper_sh, 'w') as f: f.write('#!/bin/bash\n{}\n'.format(command_line));  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)
            command_line = shell_wrapper_sh

        job = dict(target=target, executable=command_line, arguments='', queue=1, memory=memory, initial_dir=initial_dir, run_time=run_time)
        self.execute_hpc_jobs([job])


    def execute_hpc_jobs(self, hpc_jobs):
        ''' Submit given list of jobs to HPC cluster
            each job is a dictionary with following keys:
                target : arbitrary name of the job target (will be added to HPC submit file name)
                executable : path + executable name that will be run
                arguments : command line arguments for executable. May contatin unresolved template {process} that need to be filled in.
                initial_dir : path to starting dir from which command should be executed
                queue : int number of jobs to submit.
                memory : Megabytes required to be present on node for job to run,
                run_time (float) : number of hours which job allowed to run

            must accumulate cpu hours in self.cpu_usage
        '''
        must_be_implemented_in_inherited_classes
