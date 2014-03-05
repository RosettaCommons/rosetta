#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

import sys, os, re, glob, argparse, yaml
from pylab import *
from mpldatacursor import datacursor
from graphics import tango

parser = argparse.ArgumentParser()
parser.add_argument('old_benchmark')
parser.add_argument('new_benchmarks', nargs='+')
parser.add_argument('--include', '-i', default='.')
parser.add_argument('--exclude', '-x', default='^$')
parser.add_argument('--label', '-l', default='^$')
parser.add_argument('--quiet', '-q', action='store_true')
parser.add_argument('--force', '-f', action='store_true')
parser.add_argument('--output', '-o')
parser.add_argument('--dpi', '-d', type=int, default=150)
arguments = parser.parse_args()

def parse_benchmark(benchmark):
    cache_path = os.path.join(benchmark, 'percent_subangstrom.dat')

    if not arguments.force:
        try: 
            with open(cache_path) as file:
                return yaml.load(file)
        except:
            pass

    hit_percents = {}
    loop_pattern = os.path.join(benchmark, '????')
    loops = glob.glob(loop_pattern)

    for loop in loops:
        if not re.search(arguments.include, loop):
            print "Not including:", loop
            continue
        if re.search(arguments.exclude, loop):
            print "Excluding:", loop
            continue

        loop_name = os.path.basename(loop)
        hit_percents[loop_name] = parse_loop(loop)

    with open(cache_path, 'w') as file:
        string = yaml.dump(hit_percents)
        file.write(string)

    return hit_percents

def parse_loop(loop):
    job_pattern = os.path.join(loop, '*')
    jobs = glob.glob(job_pattern)

    hit_count = 0

    for job in sorted(jobs):
        rmsd = parse_job(job)
        if rmsd is None: continue
        if rmsd < 1: hit_count += 1

    print
    return hit_count / len(jobs)

def parse_job(job):
    path = os.path.join(job, 'output.txt')

    rmsd_header = 'protocols.loop_build.LoopBuildMover: loop_rms'
    score_header = 'protocols.loop_build.LoopBuildMover: total_energy'

    if not arguments.quiet:
        print 'Parsing %s...' % job
        #sys.stdout.write('\rParsing %s...' % job)
        #sys.stdout.flush()

    try:
        with open(path) as file:
            file.seek(0, 2)             # Go to the end of the file.
            file_size = file.tell()     # Figure out how big the file is.

            tail_start = max(file_size - 1000, 0)
            file.seek(tail_start, 0)
            lines = file.readlines()

        for line in lines:
            if line.startswith(rmsd_header):
                fields = line.split()
                return float(fields[-1])

    except IOError:
        print '\rNo data for:', path
        return None


def compare_benchmarks(old_benchmark, new_benchmarks, labels):
    x = y = linspace(0, 1)
    plot(x, y, '--', color=tango.grey[3])

    colors = [ tango.blue[1], tango.red[1], tango.green[2],
            tango.brown[1], tango.purple[1] ]

    for index, new_benchmark in enumerate(new_benchmarks):
        plot_benchmark(
                old_benchmark, new_benchmark, labels[index], colors[index])

    title('Percent Sub-Angstrom Decoys')
    xlabel('Standard KIC')
    ylabel('Conservative KIC')

    xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0', '0.2', '0.4', '0.6', '0.8', '1.0'])
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0', '0.2', '0.4', '0.6', '0.8', '1.0'])

    #xlim(0, 0.305)
    #ylim(0, 0.305)

    axes([0.19, 0.6, .25, .25])
    datacursor(formatter='{label}'.format)

    x = y = linspace(0, 1)
    plot(x, y, '--', color=tango.grey[3])

    for index, new_benchmark in enumerate(new_benchmarks):
        plot_benchmark(
                old_benchmark, new_benchmark, labels[index], colors[index])

    xlim(0, 0.02)
    ylim(0, 0.02)
    xticks([0, 0.02], ['0', '0.02'])
    yticks([0, 0.02], ['0', '0.02'])

    if arguments.output:
        savefig(arguments.output, dpi=arguments.dpi)
    elif not os.fork():
        show()

def plot_benchmark(old_benchmark, new_benchmark, label, color):
    old_loops = set(old_benchmark.keys())
    new_loops = set(new_benchmark.keys())

    for loop in old_loops & new_loops:
        plot(old_benchmark[loop], new_benchmark[loop], 'o',
             color=color, label='{}/{}'.format(os.path.basename(label), loop))


if __name__ == '__main__':
    old = parse_benchmark(arguments.old_benchmark)
    new = [parse_benchmark(x) for x in arguments.new_benchmarks]

    if not arguments.quiet:
        print

    if not old or not all(new):
        print "Giving up because no data was found."
        raise SystemExit

    compare_benchmarks(old, new, arguments.new_benchmarks)
