#!/usr/bin/env python

# Make sure this script isn't being run on chef.  It won't work anyway, but the 
# errors will be kinda cryptic.

from socket import gethostname
if gethostname().startswith('chef'):
    print "Don't run this script on chef."
    raise SystemExit

import os
import argparse
import subprocess
import scripting

# Specify default values for several important variables.  This section of the 
# script can be thought of as a primitive configuration file.  All of these 
# defaults can be overridden on the command line.

name_abbreviations = {
        'sandbox': 'jobs/sandbox',
        'algorithm': 'jobs/{0.algorithm}',
        'pdb': 'jobs/{0.algorithm}/{0.pdb}',
        'commit': 'jobs/{0.pdb}/{0.algorithm}/{0.commit}',
        'cluster': 'jobs/cluster/{0.algorithm}/{0.pdb}/{0.job_id:03d}',
        'small-move': 'jobs/cluster/{0.algorithm}/{0.step_size:05.2f}/{0.pdb}/{0.job_id:03d}',
        'valgrind': 'jobs/valgrind/{0.algorithm}'
}

algorithm_choices = [
        'refactor',
        'refactor-only',
        'balanced',
        'legacy',
        'loopmodel-standard',
        'loopmodel-small-move',
        'loopmodel-hybrid',
        'loopmodel-refactor',
        'loopmodel-master',
]

default_name = 'sandbox'
default_algorithm = 'refactor'
default_structure = 'structures/kic/1srp.pdb'

# Read parameters from the command line.  These options are used to setup and 
# execute a rosetta job.

class ParseIterations (argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        if values[0] == 'fast':
            values = [3, 12, 2]
            fast = True
        else:
            values = map(self.smart_int, values) + [1, 1, 1]
            fast = False

        setattr(namespace, self.dest, values[:3])
        setattr(namespace, 'fast', fast)

    @staticmethod
    def smart_int(string):
        """ Convert the given string to an integer.  Metric suffixes, if they 
        are present, will be properly translated. """
        import re

        match = re.match(r'^(\d+)([KMG]?)$', string)
        if not match:
            message = "invalid literal for smart_int() with base 10: '%s'" 
            raise ValueError(message % string)

        digits, suffix = match.groups()
        multiplier = dict(K=1e3, M=1e6, G=1e9).get(suffix, 1)

        return int(digits) * int(multiplier)


parser = argparse.ArgumentParser()

parser.add_argument('--name', '-n', default=default_name)
parser.add_argument('--description', '-m', action='store_true')
parser.add_argument('--overwrite', '-o', action='store_true')
parser.add_argument('--background', '-bg', action='store_true')
parser.add_argument('--dry-run', '-dr', action='store_true')
parser.add_argument('--debug', '-d', dest='run', action='store_false')

parser.add_argument('--structure', '-s',
        default=default_structure)
parser.add_argument('--algorithm', '-a',
        default=default_algorithm, choices=algorithm_choices)
parser.add_argument('--iterations', '-i',
        action=ParseIterations, nargs='+', default=(10, 1, 1))

parser.add_argument('--random', '-r', action='store_true')
parser.add_argument('--random-seed', '-rs', type=int)
parser.add_argument('--fewer-rotamers', '-fr', action='store_true')
parser.add_argument('--skip-centroid', '-sc', action='store_true')
parser.add_argument('--step-size', '-ss', type=float, default=20.0)
parser.add_argument('--verbose', '-v', action='store_true')
parser.add_argument('--log-mc-poses', '-lp', action='store_true')
parser.add_argument('--log-all-poses', '-lpp', action='store_true')

job = scripting.Job()
arguments = parser.parse_args()

arguments.pdb = os.path.splitext(os.path.basename(arguments.structure))[0]
arguments.commit = subprocess.check_output('git rev-parse --short HEAD', shell=True).strip()

try: arguments.job_id = int(os.environ['SGE_TASK_ID'])
except: arguments.job_id = 0

if arguments.name == 'sandbox':
    arguments.background = False
    arguments.overwrite = True

loop_path = os.path.splitext(arguments.structure)[0] + '.loop'
with open(loop_path) as file:
    arguments.loop = tuple(map(int, file.readline().split()[0:2]))

# Setup the basic job parameters, including what the job should be named, what 
# command should be executed, whether or not old jobs should be overwritten, 
# and other things like that.

def setup_hook():
    os.symlink(rosetta_subpath('apps/structures'), 'structures')
    os.symlink(rosetta_subpath('apps/kic/analysis'), 'analysis')

    if arguments.log_mc_poses or arguments.log_all_poses:
        os.mkdir('trajectory')

    from socket import gethostname
    with open('cluster.txt', 'w') as file:
        file.write(gethostname() + '\n')
        file.write(os.environ.get('JOB_ID', 'none') + '\n')

def find_rosetta_installation():
    try: 
        # This should be called from within a rosetta installation.  If this is 
        # the case, git can be used to find the root of the installation.

        with open(os.devnull, 'w') as devnull:
            from subprocess import check_output
            command = 'git', 'rev-parse', '--show-toplevel'
            output = check_output(command, stderr=devnull)
            return output.strip()

    except subprocess.CalledProcessError:
        print "No rosetta installation found!"
        raise SystemExit(1)


job.name = arguments.name
job.name.arguments = arguments
job.name.abbreviations = name_abbreviations
job.overwrite = 'yes' if arguments.overwrite else 'no'
job.foreground = not arguments.background
job.hook = setup_hook
job.dry_run = arguments.dry_run

if arguments.description:
    job.compose_description()

rosetta_path = find_rosetta_installation()
rosetta_subpath = lambda subpath: os.path.join(rosetta_path, subpath)

# Fill in the arguments to be passed into the rosetta command.  Most of the 
# arguments are given default values by this script, and do not respect the 
# default values hard-coded into rosetta.

if arguments.algorithm.startswith('loopmodel'):
    command = rosetta_subpath('source/bin/loopmodel')
    database = rosetta_subpath('database')
    loop_file = os.path.splitext(arguments.structure)[0] + '.loop'
    remodel = 'perturb_kic'
    refine = 'refine_kic'

    if arguments.algorithm == 'loopmodel-refactor':
        remodel = 'perturb_kic_refactor'
        refine = 'refine_kic_refactor'

    if arguments.skip_centroid:
        remodel = 'no'

    # This hack is not very robust, but it makes it possible to invoke the 
    # `loopmodel' executable from the master branch.  Assuming your directory 
    # structure is the same as mine.  Ha ha.

    if arguments.algorithm == 'loopmodel-master':
        command = command.replace('develop', 'master')
        database = database.replace('develop', 'master')

    job.command = command

    job.add_argument('-database %s' % database)
    job.add_argument('-in:file:s %s ' % arguments.structure)
    job.add_argument('-in:file:native %s' % arguments.structure)
    job.add_argument('-in:file:fullatom')
    job.add_argument('-out:nooutput')
    job.add_argument('-nstruct 1')
    job.add_argument('-loops:loop_file %s' % loop_file)
    job.add_argument('-loops:remodel %s' % remodel)
    job.add_argument('-loops:refine %s' % refine)
    job.add_argument('-loops:legacy_kic false')

    if arguments.algorithm == 'loopmodel-small-move':
        job.add_argument('-loops:kic_small_moves')
        job.add_argument('-loops:kic_small_move_magnitude %f' % arguments.step_size)
    if arguments.algorithm == 'loopmodel-hybrid':
        job.add_argument('-loops:kic_pivot_based')

    if arguments.iterations == [3, 12, 2]:
        job.add_argument('-loops:fast')

    job.add_argument('-mute core.io.pdb.file_data')
    job.add_argument('-unmute protocols.loop_build.LoopBuildMover')

else:
    job.command = rosetta_subpath('source/bin/KicSandbox')

    job.add_argument('-database %s' % rosetta_subpath('database'))
    job.add_argument('-kale:algorithm %s' % arguments.algorithm)
    job.add_argument('-kale:structure %s' % arguments.structure)
    job.add_argument('-kale:loop %d %d' % arguments.loop)
    job.add_argument('-kale:iterations %d %d %d' % tuple(arguments.iterations))
    job.add_argument('-kale:out:mc_poses %s' % arguments.log_mc_poses)
    job.add_argument('-kale:out:all_poses %s' % arguments.log_all_poses)

    if arguments.background:
        job.add_argument('-kale:out:progress_bar false')

if not arguments.verbose:
    job.add_argument('-mute all')
    job.add_argument('-unmute apps.pilot.kale')
    job.add_argument('-unmute protocols.loop_build.LoopBuildMover')

if not arguments.fewer_rotamers:
    job.add_argument('-ex1 -ex2 -extrachi_cutoff 0')

if arguments.random_seed:
    job.add_argument('-constant_seed')
    job.add_argument('-jran %d' % arguments.random_seed)

elif not arguments.random:
    job.add_argument('-constant_seed')

# Either run or debug the job, depending on what was requested.  

if arguments.run: job.run()
else: job.debug()

