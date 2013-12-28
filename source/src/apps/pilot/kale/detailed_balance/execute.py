#!/usr/bin/env python

import os
import argparse
import scripting

# Specify default values for several important variables.  This section of the 
# script can be thought of as a primitive configuration file.  All of these 
# defaults can be overridden on the command line.

name_abbreviations = {
        'sandbox':       'jobs/sandbox',
}

default_name = 'sandbox'
default_structure = 'structures/contrived/valine.6.pdb'
default_residues = {
        'structures/contrived/valine.6.pdb': (3, 3),
        'structures/contrived/lysine.6.pdb': (3, 3),
        'structures/benchmark/1srp.pdb': (308, 319)
}

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

parser.add_argument('--structure', '-s', default=default_structure)
parser.add_argument('--residues', '-x', type=int, nargs=2, default=None)
parser.add_argument('--iterations', '-i', type=smart_int, default=10000)
parser.add_argument('--random', '-r', action='store_true')
parser.add_argument('--random-seed', '-rs', type=int)
parser.add_argument('--verbose', '-v', action='store_true')
parser.add_argument('--dump-pose', '-dp', type=int, default=0)
parser.add_argument('--no-progress-bar', '-npb', action='store_true')

job = scripting.Job()
arguments = parser.parse_args()

if arguments.residues is None:
    try: arguments.residues = default_residues[arguments.structure]
    except KeyError:
        print "Residues not known for structure '%s'." % arguments.structure
        raise SystemExit

arguments.pdb_string = os.path.basename(arguments.structure)

# Setup the basic job parameters, including what the job should be named, what 
# command should be executed, whether or not old jobs should be overwritten, 
# and other things like that.

# rosetta_subpath (fold)
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

job.command = rosetta_subpath('source/bin/sidechain_tests')
job.overwrite = 'yes' if arguments.overwrite else 'no'
job.foreground = not arguments.background
job.hook = setup_hook
job.dry_run = arguments.dry_run

if arguments.description:
    job.compose_description()

# Fill in the arguments to be passed into the sampling command.  Most of the 
# arguments are given default values by this script, and do not respect the 
# default values hard-coded into rosetta.

job.add_argument('-kale:in:pdb %s' % arguments.structure)
job.add_argument('-kale:in:residues %d %d' % tuple(arguments.residues))
job.add_argument('-kale:in:iterations %d' % arguments.iterations)
job.add_argument('-kale:out:trajectory %s' % arguments.dump_pose)

if arguments.no_progress_bar or arguments.background:
    job.add_argument('-kale:out:quiet true')

if not arguments.verbose:
    job.add_arguments('-mute all', '-unmute apps.pilot.kale')

if arguments.random_seed:
    job.add_argument('-constant_seed')
    job.add_argument('-jran %d' % arguments.random_seed)

elif not arguments.random:
    job.add_argument('-constant_seed')

# Either run or debug the job, depending on what was requested.  

if arguments.run:
    os.nice(10)
    job.run()
else:
    job.debug()

