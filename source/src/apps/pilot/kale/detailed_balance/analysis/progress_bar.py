#!/usr/bin/env python
# encoding: utf-8

import os, sys
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('job')
arguments = parser.parse_args()

def main():
    # Figure out how many iterations there will be.
    job = os.path.join(arguments.job, 'job.sh')
    with open(job) as file:
        for line in file:
            fields = line.strip().split()
            if fields and fields[0] == '-kale:mc:iterations':
                total_iterations = int(fields[1])

    # Figure out how many iterations have already happened.
    coordinates = os.path.join(arguments.job, 'coordinates.dat')

    while True:
        time.sleep(0.02)
        completed_iterations = os.path.getsize(coordinates) // 182
        progress = completed_iterations, total_iterations

        sys.stdout.write('\033[?25l\r[{0}/{1}]'.format(*progress))
        sys.stdout.flush()

        # Stop once the job finished.
        if completed_iterations >= total_iterations:
            sys.stdout.write('\033[?25l\r[{1}/{1}]\n'.format(*progress))
            break

if __name__ == '__main__':
    try: main()
    except KeyboardInterrupt: print
    finally: sys.stdout.write('\033[?25h')

