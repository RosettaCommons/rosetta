# -*- coding: utf-8 -*-
# :noTabs=true:

import time as time_module

import sys, signal, imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) + '/base.py')  # A bit of Python magic here, what we trying to say is this: from base import *, but path to base is calculated from our source location  # from base import HPC_Driver, execute, NT


class MultiCore_HPC_Driver(HPC_Driver):
    # def execute(self, command_line, working_dir, target, memory=256, run_time=24):
    #     self.cpu_usage -= time.time()/60./60.
    #     log = execute('Executing {}'.format(target), 'cd {} && {}'.format(working_dir, command_line), tracer=self.tracer, return_='output')
    #     with file(self.working_dir+'/.hpc.{target}.log'.format(target=target), 'w') as f: f.write(log)
    #     self.cpu_usage += time.time()/60./60.


    def submit_hpc_job(self, name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory=512, time=12, block=True):
        cpu_usage = -time_module.time()/60./60.

        _jobs_ = []
        def mfork():
            ''' Check if number of child process is below cpu_count. And if it is - fork the new pocees and return its pid.
            '''
            while len(_jobs_) >= self.cpu_count:
                for p in _jobs_[:] :
                    r = os.waitpid(p, os.WNOHANG)
                    if r == (p, 0):  # process have ended without error
                        _jobs_.remove(p)
                    elif r[0] == p :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                        for p in _jobs_: os.waitpid(p, 0)
                        self.trace('Some of the HPC jobs terminated abnormally!')
                        sys.exit(1)

                if len(_jobs_) >= jobs: time.sleep(.5)

            sys.stdout.flush()
            pid = os.fork()
            if pid: _jobs_.append(pid) # We are parent!
            return pid


        previous_handler = signal.getsignal(signal.SIGINT)

        def signal_handler(signal_, f):
            self.trace('Ctrl-C pressed... killing child jobs...')
            for j in _jobs_:
                os.killpg(os.getpgid(j), signal.SIGKILL)

            if type(previous_handler) == function: previous_handler(signal_, f)

        signal.signal(signal.SIGINT, signal_handler)

        process = 0
        for i in range(jobs_to_queue):

            pid = mfork()
            if not pid:  # we are child process
                command_line = 'cd {} && {} {}'.format(working_dir, executable, arguments.format(process=process) )
                log = execute('Running job {}.{}...'.format(name, i), command_line, tracer=self.tracer, return_='output')
                with file(self.working_dir+'/hpc.{name}.{i}.log'.format(**vars()), 'w') as f: f.write(command_line+'\n'+log)
                sys.exit(0)

            process += 1

        for p in _jobs_: os.waitpid(p, 0)  # waiting for all child process to termintate...

        cpu_usage += time_module.time()/60./60.

        self.cpu_usage += cpu_usage * jobs_to_queue  # approximation...

        signal.signal(signal.SIGINT, previous_handler)
