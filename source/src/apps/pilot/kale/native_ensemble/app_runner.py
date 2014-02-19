#!/usr/bin/env python

import argparse, copy
import os, sys, subprocess, platform

parser = argparse.ArgumentParser()
parser.add_argument('--build', '-b', default='release')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--algorithm', '-a', default='walking')
parser.add_argument('--iterations', '-i', type=int, default=500)
parser.add_argument('--temp-levels', '-j', type=int, default=4)
parser.add_argument('--temp-max', '-tx', type=float, default=2)
parser.add_argument('--temp-freq', '-tf', type=int, default=1000)
parser.add_argument('--weights', '-w', nargs=2, default=(1, 1))
parser.add_argument('--frequency', '-f', type=int, default=10)
parser.add_argument('--message', '-m', default='')
parser.add_argument('--truncate', '-0', action='store_true')
parser.add_argument('--runner', '-r', default='simple')
parser.add_argument('--database', '-d', default='sandbox.db')
arguments = parser.parse_args()

options = copy.deepcopy(arguments)
options.weights = '{} {}'.format(*arguments.weights)
options.build_flags = ''
options.execute_flags = '-d' if arguments.debug else ''

if arguments.truncate:
    print 'Truncating %s' % arguments.database
    os.remove(arguments.database)

if arguments.runner == 'mysql':     # (fold)
    if platform.node() == 'glycine':
        print "Must run in sandbox mode on glycine."
        raise SystemExit

    command = '''\
    rosetta_build --build=mpi && \\
    rosetta_execute native_ensemble -v -- \\
        -out:use_database \\
        -inout:dbms:mode mysql \\
        -inout:dbms:database_name kale \\
        -inout:dbms:user kale \\
        -inout:dbms:port 3306 \\
        -inout:dbms:password $(abraxas mysql) \\
        -inout:dbms:host guybrush.ucsf.edu \\
        -s structures/benchmark/154l.pdb \\
        -native_ensemble:loop 153 164 \\
        -native_ensemble:weights {0.weights} \\
        -native_ensemble:algorithm {0.algorithm} \\
        -native_ensemble:iterations {0.iterations} \\
        -native_ensemble:frequency {0.frequency} \\
        -native_ensemble:message "'{0.message}'"
    '''.format(options)

elif arguments.runner == 'mpi':     # (fold)
    command = '''\
    rosetta_build --build=mpi && \\
    mpirun -np {0.temp_levels} \\
      ~/rosetta/develop/source/bin/native_ensemble \\
        -no_output \\
        -database_name {0.database} \\
        -s structures/benchmark/154l.pdb \\
        -native_ensemble:loop 153 164 \\
        -native_ensemble:weights {0.weights} \\
        -native_ensemble:algorithm {0.algorithm} \\
        -native_ensemble:iterations {0.iterations} \\
        -native_ensemble:frequency {0.frequency} \\
        -native_ensemble:message "'{0.message}'" \\
        -tempering:temp:levels {0.temp_levels}
    '''.format(options)

elif arguments.runner == 'guybrush':     # (fold)
    #-inout:dbms:password $(abraxas mysql) \\
    command = '''\
    mpirun -np {0.temp_levels} \\
      ~/rosetta/develop/source/bin/native_ensemble.mpimysql.linuxgccrelease \\
        -no_output \\
        -inout:dbms:mode mysql \\
        -inout:dbms:database_name kale \\
        -inout:dbms:user kale \\
        -inout:dbms:port 3306 \\
        -inout:dbms:host guybrush.ucsf.edu \\
        -inout:dbms:password pJarth26xSRs2RrC \\
        -s structures/benchmark/154l.pdb \\
        -native_ensemble:loop 153 164 \\
        -native_ensemble:weights {0.weights} \\
        -native_ensemble:algorithm walking/para-temp \\
        -native_ensemble:iterations {0.iterations} \\
        -native_ensemble:frequency {0.frequency} \\
        -native_ensemble:message "'{0.message}'" \\
        -tempering:temp:low 1.0 \\
        -tempering:temp:high {0.temp_max} \\
        -tempering:temp:levels {0.temp_levels} \\
        -tempering:stride {0.temp_freq}
    '''.format(options)

elif arguments.runner == 'callgrind':  # (fold)
    command = '''\
    rosetta_build --build={0.build} {0.build_flags} && \\
    valgrind --tool=callgrind \\
      ~/rosetta/develop/source/bin/native_ensemble \\
        -database_name {0.database} \\
        -s structures/benchmark/154l.pdb \\
        -native_ensemble:loop 153 164 \\
        -native_ensemble:weights {0.weights} \\
        -native_ensemble:algorithm {0.algorithm} \\
        -native_ensemble:iterations {0.iterations} \\
        -native_ensemble:frequency {0.frequency} \\
        -native_ensemble:message "'{0.message}'" \\
        -tempering:temp:low {0.temp_min} \\
        -tempering:temp:high {0.temp_max} \\
        -tempering:temp:levels {0.temp_levels} \\
        -tempering:stride {0.temp_stride}
    '''.format(options)

elif arguments.runner == 'simple':  # (fold)
    command = '''\
    rosetta_build --build={0.build} {0.build_flags} && \\
    rosetta_execute -v {0.execute_flags} -- \\
        -database_name {0.database} \\
        -s structures/benchmark/154l.pdb \\
        -native_ensemble:loop 153 164 \\
        -native_ensemble:weights {0.weights} \\
        -native_ensemble:algorithm {0.algorithm} \\
        -native_ensemble:iterations {0.iterations} \\
        -native_ensemble:frequency {0.frequency} \\
        -native_ensemble:message "'{0.message}'" \\
        -tempering:temp:low {0.temp_min} \\
        -tempering:temp:high {0.temp_max} \\
        -tempering:temp:levels {0.temp_levels} \\
        -tempering:stride {0.temp_stride}
    '''.format(options)

else:   # (fold)
    print "Unknown runner '%s'." % arguments.runner
    sys.exit()


subprocess.call(command, shell=True)
