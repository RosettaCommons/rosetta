#!/bin/env python
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   performance-benchmark.py
## @brief  Script for running mini benchmarks
## @author Sergey Lyskov


import os, commands, sys, re, subprocess


# Try to identity plaform by using scons compiliation feature.
def getPlatformID():
    print "Identifying platform...\n"
    pl = commands.getoutput("""scons Abracadabra log=platform mode=release""")
    lines = pl.split('\n')
    for s in lines:
        if  len( s.split() ) > 1 and s.split()[0] == 'Platform:':
            platform = s.split()[1]
            print "Platform found: " + platform
            return platform
    return "PlatformWasNotFound!!!"  # <-- That should not reall happend.


def main(args):
    database="-database ~/minirosetta_database"

    # The flags that you must have in order for the performance benchmark to work.
    # Here we currently have the flags required for the ligand docking benchmark to work.
    required_flags = "-in:file:extra_res_path extra_params"

    cl=required_flags
    if len(args) <= 1:
        print "No database path supplied... usign defaut one."
    else:
        database = " ".join(args[1:])
        cl += " -mute core protocols " + database
    print "Comand line arguments:", cl

    #platform = getPlatformID()
    #f = os.popen("""cd demo && ./../build/demo/%s/benchmark %s 1>&2""" % (platform, cl), 'r')
    f = os.popen("""cd src/apps/benchmark/performance && ./../../../../bin/performance_benchmark.default.linuxgccrelease %s 1>&2""" % (cl), 'r')

    #f = subprocess.Popen("""cd demo && ./../build/demo/%s/benchmark %s""" % (platform, cl), bufsize=0, shell=True, stdout=subprocess.PIPE).stdout
    for line in f:
        print line,
        #sys.stdout.flush()
    f.close()

    print "Done!"



if __name__ == "__main__": main(sys.argv)

