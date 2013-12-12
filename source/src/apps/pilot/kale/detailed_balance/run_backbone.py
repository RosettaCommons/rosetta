#!/usr/bin/env python

import os
import argparse
import scripting

# Specify default values for several important variables.  This section of the 
# script can be thought of as a primitive configuration file.  All of these 
# defaults can be overridden on the command line.

name_abbreviations = {
        'sm':            'jobs/sampling/{0.sampling_move}',
        'sampling':      'jobs/sampling/{0.sampling_move}',
        'free':          'jobs/{0.sampling_move}-{0.score_function}/free/{0.closure_move}',
        'fixed':         'jobs/{0.sampling_move}-{0.score_function}/fixed/{0.closure_move}',
        'fixed-pivots':  'jobs/{0.sampling_move}-{0.score_function}/fixed/{0.closure_move}/{0.pivot_string}',
        'sandbox':       'jobs/sandbox',
}
names = name_abbreviations.keys()
default_name = 'sandbox'

closure_moves = 'naive', 'balanced',
sampling_moves = 'uniform', 'rama', 'walking', 'vicinity'
breadth_moves = 'uniform', 'omega', 'rama'
score_functions = 'rama', 'off'

default_closure_move = 'naive'
default_sampling_move = 'uniform'
default_breadth_move = 'omega'
default_score_function = 'off'

default_structure = 'structures/linear/6.symmetry.pdb'
default_pivots = {
        'structures/linear/4.helix.pdb': (1, 2, 4),
        'structures/linear/5.helix.pdb': (2, 3, 4),
        'structures/linear/6.symmetry.pdb': (2, 3, 5),
        'structures/linear/7.helix.pdb': (2, 4, 6),
        'structures/linear/7.1ubq.pdb': (2, 4, 6),
        'structures/linear/14.1eco.pdb': (2, 7, 12),
        'structures/linear/14.cycle.pdb': (2, 7, 12),
        'structures/kic/1srp.pdb': (308, 319)
}

try: 
    with open(os.devnull, 'w') as devnull:
        from subprocess import check_output
        command = 'git', 'rev-parse', '--show-toplevel'
        output = check_output(command, stderr=devnull)
        rosetta_path = output.strip()

except subprocess.CalledProcessError:
    print "No rosetta installation found!"
    raise SystemExit(1)

rosetta_subpath = lambda subpath: os.path.join(rosetta_path, subpath)

def smart_int(string):
    """ Convert the given string to an integer.  Metric suffixes, if they are 
    present, will be properly translated. """
    import re

    match = re.match(r'^(\d+)([KMG]?)$', string)
    if not match:
        message = "invalid literal for smart_int() with base 10: '%s'" % string
        raise ValueError(message)

    digits, suffix = match.groups()
    multiplier = dict(K=1e3, M=1e6, G=1e9).get(suffix, 1)

    return int(digits) * int(multiplier)


# Read parameters from the command line.  These options are used to setup and 
# execute a rosetta job.

parser = argparse.ArgumentParser()

parser.add_argument('--name', '-n', default=default_name)
parser.add_argument('--description', '-m', action='store_true')
parser.add_argument('--overwrite', '-o', action='store_true')
parser.add_argument('--background', '-bg', action='store_true')
parser.add_argument('--dry-run', '-z', action='store_true')
parser.add_argument('--debug', '-d', dest='run', action='store_false')

parser.add_argument('--closure-move', '-cm',
        default=default_closure_move, choices=closure_moves)
parser.add_argument('--sampling-move', '-sm',
        default=default_sampling_move, choices=sampling_moves)
parser.add_argument('--breadth-move', '-bm',
        default=default_breadth_move, choices=breadth_moves)
parser.add_argument('--score-function', '-sf',
        default=default_score_function, choices=score_functions)

parser.add_argument('--structure', '-s', default=default_structure)
parser.add_argument('--pivots', '-p', type=int, nargs=3, default=None)
parser.add_argument('--iterations', '-i', type=smart_int, default=1000)
parser.add_argument('--temperature', '-t', type=float, default=1)
parser.add_argument('--weights', '-w', type=float, nargs=3, default=(1,0,0))
parser.add_argument('--random', '-r', action='store_true')
parser.add_argument('--random-seed', '-rs', type=int)
parser.add_argument('--verbose', '-v', action='store_true')
parser.add_argument('--dump-pose', '-dp', type=int, default=0)
parser.add_argument('--dump-stats', '-ds', type=int, default=1)
parser.add_argument('--no-progress-bar', '-npb', action='store_true')

job = scripting.Job()
arguments = parser.parse_args()

# Setup the basic job parameters, including what the job should be named, what 
# command should be executed, whether or not old jobs should be overwritten, 
# and other things like that.

def setup_hook():
    os.mkdir('trajectory')
    os.symlink(rosetta_subpath('apps/structures'), 'structures')
    os.symlink(rosetta_subpath('apps/monte_carlo/analysis'), 'analysis')


job.name = arguments.name
job.name.arguments = arguments
job.name.abbreviations = name_abbreviations

if arguments.name == 'sandbox':
    arguments.background = False
    arguments.overwrite = True

job.command = rosetta_subpath('source/bin/sampling_test')
job.overwrite = 'yes' if arguments.overwrite else 'no'
job.foreground = not arguments.background
job.hook = setup_hook
job.dry_run = arguments.dry_run

if arguments.description:
    job.compose_description()

# Fill in the arguments to be passed into the sampling command.  Most of the 
# arguments are given default values by this script, and do not respect the 
# default values hard-coded into rosetta.

if arguments.pivots is None:
    try: arguments.pivots = default_pivots[arguments.structure]
    except KeyError:
        print "Pivots not known for structure '%s'." % arguments.structure
        raise SystemExit

job.add_argument('-kale:in:pdb %s' % arguments.structure)
job.add_argument('-kale:kic:pivots %d %d %d' % tuple(arguments.pivots))
job.add_argument('-kale:kic:closure_move %s' % arguments.closure_move)
job.add_argument('-kale:kic:sampling_move %s' % arguments.sampling_move)
job.add_argument('-kale:kic:breadth_move %s' % arguments.breadth_move)
job.add_argument('-kale:kic:move_weights %f %f %f' % tuple(arguments.weights))
job.add_argument('-kale:kic:score_function %s' % arguments.score_function)
job.add_argument('-kale:kic:dump_pose %s' % arguments.dump_pose)
job.add_argument('-kale:kic:dump_stats %s' % arguments.dump_stats)
job.add_argument('-kale:mc:iterations %d' % arguments.iterations)
job.add_argument('-kale:mc:temperature %d' % arguments.temperature)

if arguments.random_seed:
    job.add_argument('-constant_seed')
    job.add_argument('-jran %d' % arguments.random_seed)

elif not arguments.random:
    job.add_argument('-constant_seed')

if not arguments.verbose:
    job.add_arguments('-mute all', '-unmute apps.pilot.kale')

if arguments.no_progress_bar or arguments.background:
    job.add_argument('-kale:out:quiet true')

arguments.pdb_string = os.path.basename(arguments.structure)
arguments.pivot_string = '%d-%d-%d' % tuple(arguments.pivots)

# Either run or debug the job, depending on what was requested.  

if arguments.run:
    os.nice(10)
    job.run()
else:
    job.debug()

