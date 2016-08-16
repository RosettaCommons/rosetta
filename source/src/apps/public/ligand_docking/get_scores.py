#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# Suggested Usage:
# for f in *.out; do echo $f; ../get_scores.py < $f > ${f:0:4}_suffix.tab; done

import sys
from optparse import OptionParser


def main(argv):
    '''
    Converts the SCORES lines from an atomtree diff silent output file
    into a tabular format that can be read into R, Excel, etc.
    '''
    parser = OptionParser(usage="usage: %prog [SILENT.OUT]")
    parser.set_description(main.__doc__)
    # parser.add_option("-short", ["--long"],
    #   action="store|store_true|store_false",
    #   default=True|False|...
    #   type="string|int|float",
    #   dest="opt_name",
    #   help="store value in PLACE",
    #   metavar="PLACE",
    # )
    (options, args) = parser.parse_args(args=argv)

    if len(args) == 0:
        infile = sys.stdin
        outfile = sys.stdout
    elif len(args) == 1:
        infile = open(args[0])
        outfile = sys.stdout
    else:
        parser.print_help()
        print "Too many arguments!"
        return 1

    first = True
    for line in infile:
        if not line.startswith("SCORES "): continue
        if "is_reference_pose 1" in line: continue # skip ref. structures
        f = line.rstrip().split()
        if first:
            first = False
            i = 0
            while i < len(f):
                print f[i],
                i += 2
            print
        i = 1
        while i < len(f):
            print f[i],
            i += 2
        print

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
