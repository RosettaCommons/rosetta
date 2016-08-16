#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# Suggested usage:
# for f in *silent.out; do ~/rosetta/rosetta_source/bin/extract_atomtree_diffs.linuxiccrelease @EXTRACT.txt -s $f -tags $(../best_ifaceE.py -n3 $f); done

import sys
from optparse import OptionParser

def main(argv):
    """
    Prints the tag(s) of the "best" structures in the
    provided silent output file (atomtree diff format).
    By default, it takes the top 5% of structures by total score,
    sorts those by the energy across the binding interface,
    and prints out the top 10.

    In the Bash shell, this list of tags can be used directly
    in the extract_atomtree_diffs command line:

        extract_atomtree_diffs ... -s silent.out -tags $(best_ifaceE.py silent.out)

    """
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
    parser.add_option("-n",
      default=10,
      type="int",
      help="print tags for top N structures",
      metavar="N",
    )
    parser.add_option("-f", "--frac",
      default=0.05,
      type="float",
      help="keep only top X structures by total score",
      metavar="X",
    )
    parser.add_option("-t", "--total-score",
      default=False,
      action='store_true',
      help="sort by total score instead of interface delta",
    )
    parser.add_option("-p", "--print-score",
      default=False,
      action='store_true',
      help="print scores alongside structure tags",
    )
    parser.add_option("--no-fa-pair",
      default=False,
      action='store_true',
      help="exclude fa_pair from interface_delta, as for Davis & Baker 2008 JMB",
    )
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

    # Extract score info from file
    colnames = {} # str -> int
    structs = []
    first = True
    has_atom_pair_constr = False
    has_fa_pair = False
    for line in infile:
        if not line.startswith("SCORES "): continue
        if "is_reference_pose 1" in line: continue # skip ref. structures
        f = line.rstrip().split()
        # Set up dictionary of column names
        if first:
            first = False
            f[0] = "tag" # change "SCORES" to "tag"
            for i in xrange(0,len(f),2):
                colnames[f[i]] = i+1
            has_atom_pair_constr = "if_atom_pair_constraint" in colnames
            has_fa_pair = "if_fa_pair" in colnames and options.no_fa_pair
        # Select desired scores
        ligand_is_touching = int(f[colnames["ligand_is_touching"]])
        if options.total_score:
            score = float(f[colnames["total_score"]])
        else:
            score = float(f[colnames["interface_delta"]])
            if has_atom_pair_constr: score -= float(f[colnames["if_atom_pair_constraint"]])
            if has_fa_pair: score -= float(f[colnames["if_fa_pair"]])
        if ligand_is_touching:
            structs.append((f[colnames["tag"]], float(f[colnames["total_score"]]), score))

    # Sort and select
    structs.sort(lambda a,b: cmp(a[1],b[1])) # total score, smallest (most neg.) values first
    structs = structs[:int(round(options.frac*len(structs)))] # top 5%
    structs.sort(lambda a,b: cmp(a[2],b[2])) # interface energy, smallest (most neg.) values first
    for s in structs[:min(len(structs), options.n)]:
        if options.print_score: print s[0], s[2]
        else: print s[0] # tag
        #print s # debugging

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
