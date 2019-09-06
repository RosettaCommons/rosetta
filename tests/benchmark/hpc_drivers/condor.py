# -*- coding: utf-8 -*-
# :noTabs=true:

import os, sys, time
import stat as stat_module

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) + '/base.py')  # A bit of Python magic here, what we trying to say is this: from base import *, but path to base is calculated from our source location  # from base import HPC_Driver, execute, NT


T_condor_job_template = '''
universe     = {universe}
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

# node requirelement
Requirements = {requirements}

machine_count = {machine_count}
queue {jobs_to_queue}
'''
#universe     = vanilla
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
    def head_node_execute(self, message, command_line, *args, **kwargs):
        head_node = self.config['condor']['head_node']

        command_line = f"ssh {head_node} cd `pwd` '&& {command_line}'"
        return execute(f'Executiong on {head_node}: {message}' if message else '', command_line, *args, **kwargs)


    def get_condor_accumulated_usage(self, user='$USER@'):
        # Expected output: sergey@UT64 0.50 0.50 1.00 0 196.86 12/26/2010 23:30  9/27/2011 23:55
        for _ in range(8):
            try:
                o = self.head_node_execute('', 'condor_userprio -all -allusers | grep {}'.format(user), return_='output', terminate_on_failure=False, silent=True).split()
                if len(o) >= 6: return max(1.0, float( o[5] ) )
                else: return 0.0
            except ValueError as _: pass

            time.sleep(32)

        return 0.0


    @property
    def number_of_cpu_per_node(self): return int( self.config['condor']['mpi_cpu_per_node'] )

    @property
    def maximum_number_of_mpi_cpu(self):
        return self.number_of_cpu_per_node * int( self.config['condor']['mpi_maximum_number_of_nodes'] )

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
        self.head_node_execute('Releasing condor jobs...', 'condor_release $USER', return_='tuple', silent=True)

        s = self.head_node_execute('', f'condor_q $USER | grep $USER | grep {condor_job_id}', return_='output', terminate_on_failure=False, silent=True)
        if s: return False
        else:
            #self.tracer('Waiting for condor to finish the jobs... DONE')
            self.jobs.remove(condor_job_id)
            self.cpu_usage += self.get_condor_accumulated_usage()
            return True  # jobs already finished, we return empty list to prevent double counting of cpu_usage


    def cancel_job(self, condor_job_id):
        self.head_node_execute(f'Condor_HPC_Driver.canceling job {condor_job_id}...', f'condor_rm {condor_job_id}', terminate_on_failure=False)


    def submit_hpc_job(self, name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory=512, time=12, block=True, shell_wrapper=False):
        print('submit_hpc_job is DEPRECATED and will be removed in near future, please use submit_serial_hpc_job  instead!')
        return self.submit_serial_hpc_job(name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory, time, block, shell_wrapper)


    def submit_serial_hpc_job(self, name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory=512, time=12, block=True, shell_wrapper=False):
        self.cpu_usage -= self.get_condor_accumulated_usage()

        # creating shell wrapper in order to reliably capture output
        #execute_sh = os.path.abspath( self.working_dir + '/.hpc.execute.{}.sh'.format(name) )
        #with open(execute_sh, 'w') as f: f.write('#!/bin/bash\n$*\n');  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)

        arguments = arguments.format(process='$(Process)')
        run_time=int(time*60*60)

        if shell_wrapper:
            shell_wrapper_sh = os.path.abspath(self.working_dir + f'/hpc.{name}.shell_wrapper.sh')
            with open(shell_wrapper_sh, 'w') as f: f.write('#!/bin/bash\n{} {}\n'.format(executable, arguments));  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)
            executable, arguments = shell_wrapper_sh, ''

        #jobs = []
        #???? p = dict(j, , arguments=j['arguments'].format(process='$(Process)') )

        #execute_sh = os.path.abspath(self.working_dir + '/.hpc.execute.{}.sh'.format(j['target']))
        #with open(execute_sh, 'w') as f: f.write('#!/bin/bash\n{executable} {arguments}\n'.format(**j));  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)

        condor_file = working_dir + '/.hpc.{}.condor'.format(name)
        condor_spec = T_condor_job_template.format(universe='vanilla', name=name, executable=executable,
                                                   arguments=arguments,
                                                   working_dir=working_dir, log_dir=log_dir,
                                                   memory=memory, process='$(Process)', run_time=run_time,

                                                   jobs_to_queue=jobs_to_queue,
                                                   machine_count = 1,
                                                   requirements = '',
        )
                                                   #requirements=self.config.get('condor', 'requirements', ''),
                                                   #execute_sh=execute_sh)

        with open(condor_file, 'w') as f: f.write(condor_spec)
        condor_job_id = int( self.head_node_execute('Submitting jobs to condor...', 'cd {} && condor_submit {}'.format(self.working_dir, condor_file),
                                     tracer=self.tracer, return_='output').split()[-1][:-1] )

          #jobs.append( NT(condor_file=condor_file, condor_job_id=condor_job_id, ) )

        self.jobs.append(condor_job_id)

        if block:
            self.wait_until_complete( [condor_job_id] )
            return None

        else: return condor_job_id





    def submit_mpi_hpc_job(self, name, executable, arguments, working_dir, log_dir, memory=512, time=12, block=True, process_coefficient="1", requested_nodes=1, requested_processes_per_node=1):
        ''' submit jobs as MPI job
            process_coefficient should be string representing fraction of process to launch on each node, for example '3 / 4' will start only 75% of MPI process's on each node
        '''
        max_cpu = self.maximum_number_of_mpi_cpu
        if requested_nodes*requested_processes_per_node > max_cpu : raise Exception(f'Condor_HPC_Driver.submit_mpi_hpc_job: requested_nodes times requested_processes_per_node should be below `maximum_number_of_mpi_cpu` (got requested_nodes={requested_nodes}, requested_processes_per_node={requested_processes_per_node} when maximum_number_of_mpi_cpu={self.maximum_number_of_mpi_cpu()})')

        self.cpu_usage -= self.get_condor_accumulated_usage()

        # creating shell wrapper in order to reliably capture output
        #execute_sh = os.path.abspath( self.working_dir + '/.hpc.execute.{}.sh'.format(name) )
        #with open(execute_sh, 'w') as f: f.write('#!/bin/bash\n$*\n');  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)

        arguments = arguments.format(process='$(Process)')
        run_time=int(time*60*60)

        with open( os.path.dirname(__file__) + '/condor_open_mpi_script.template.sh' ) as f: mpi_script_template = f.read()

        mpi_script = f'{working_dir}/.hpc.condor_open_mpi_script.sh'


        with open(mpi_script, 'w') as f: f.write( mpi_script_template.format(process_coefficient=process_coefficient) )
        os.chmod(mpi_script, stat_module.S_IRUSR | stat_module.S_IWUSR | stat_module.S_IXUSR)

        #jobs = []
        #???? p = dict(j, , arguments=j['arguments'].format(process='$(Process)') )

        #execute_sh = os.path.abspath(self.working_dir + '/.hpc.execute.{}.sh'.format(j['target']))
        #with open(execute_sh, 'w') as f: f.write('#!/bin/bash\n{executable} {arguments}\n'.format(**j));  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)

        condor_file = working_dir + '/.hpc.{}.condor'.format(name)

        condor_spec = T_condor_job_template.format(universe='parallel', name=name,
                                                   executable=mpi_script, arguments = executable + ' ' + arguments,
                                                   working_dir=working_dir, log_dir=log_dir,
                                                   memory=memory, process='$(Process)', run_time=run_time,
                                                   jobs_to_queue = 1,
                                                   machine_count = requested_nodes*requested_processes_per_node,
                                                   requirements = self.config['condor']['mpi_requirements'],
        )

        with open(condor_file, 'w') as f: f.write(condor_spec)
        condor_job_id = int( self.head_node_execute('Submitting jobs to condor...', 'cd {} && condor_submit {}'.format(self.working_dir, condor_file), tracer=self.tracer, return_='output').split()[-1][:-1] )

        #jobs.append( NT(condor_file=condor_file, condor_job_id=condor_job_id, ) )

        self.jobs.append(condor_job_id)

        if block:
            self.wait_until_complete( [condor_job_id] )
            return None

        else: return condor_job_id
