#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import sys
from optparse import OptionParser

def pair_up(itr):
    '''Generates (key,val) pairs from [key1, val1, key2, val2, ...].'''
    i = iter(itr)
    while True: yield (i.next(), i.next())

def main(argv):
    '''
    Parse an "atomtree diff" format silent file and remove all but the top 5% by total score.
    Also remove duplicate reference structures, which makes them faster to load
    in Mini because fewer reference Pose objects have to be held in memory.
    '''
    parser = OptionParser(usage="usage: %prog IN > OUT")
    parser.set_description(main.__doc__)
    # parser.add_option("-short", ["--long"],
    #   action="store|store_true|store_false",
    #   default=True|False|...
    #   type="string|int|float",
    #   dest="opt_name",
    #   help="store value in PLACE",
    #   metavar="PLACE",
    # )
    parser.add_option("-f", "--frac",
      default=0.05,
      type="float",
      help="fraction of poses to retain (default 0.05)",
      metavar="X",
    )
    (options, args) = parser.parse_args(args=argv)

    if len(args) == 1:
        infile = args[0]
    else:
        parser.print_help()
        print "Too many/few arguments!"
        return 1

    # Extract score info from file
    colnames = {} # str -> int
    structs = []
    first = True
    has_atom_pair_constr = False
    for line in open(infile):
        if not line.startswith("SCORES "): continue
        if "is_reference_pose 1" in line: continue # skip ref. structures
        f = line.rstrip().split()
        # Set up dictionary of column names
        if first:
            first = False
            f[0] = "tag" # change "SCORES" to "tag"
            for i in xrange(0,len(f),2):
                colnames[f[i]] = i+1
        # Select desired scores
        ligand_is_touching = int(f[colnames["ligand_is_touching"]])
        score = float(f[colnames["total_score"]])
        if ligand_is_touching:
            structs.append((f[colnames["tag"]], score))

    # Sort and select
    structs.sort(lambda a,b: cmp(a[1],b[1])) # total score, smallest (most neg.) values first
    structs = structs[:int(round(options.frac*len(structs)))] # top 5%
    struct_set = set([s[0] for s in structs])
    assert len(structs) == len(struct_set)

    curr_tag = None
    curr_scores = None
    curr_block = None
    ref_block = None
    ref_tag = None

    # Process one POSE_TAG / END_POSE_TAG block at a time.
    # If we encounter a block that looks just like the most
    # recent reference pose, then we discard it.
    out = sys.stdout
    for line in open(infile):
        if line.startswith("POSE_TAG"):
            f = line.split()
            if curr_tag is not None:
                raise ValueError("Can't start new pose %s without ending %s" % (f[1], curr_tag))
            curr_tag = f[1]
            curr_scores = None
            curr_block = []
            curr_block.append(line)
        elif line.startswith("SCORES"):
            curr_scores = dict(pair_up(line.split()))
            curr_block.append(line)
        elif line.startswith("END_POSE_TAG"):
            f = line.split()
            if f[1] != curr_tag:
                raise ValueError("Can't end pose %s without ending %s first" % (f[1], curr_tag))
            curr_block.append(line)
            write_me = True
            if 'is_reference_pose' in curr_scores and curr_scores['is_reference_pose'] == '1':
                if ref_block is not None:
                    # Have to change tag names in order for the blocks to compare equal
                    ref_block = [l.replace(ref_tag, curr_tag) for l in ref_block]
                    write_me = (curr_block != ref_block)
                ref_tag = curr_tag
                ref_block = curr_block
            else:
                write_me = (curr_tag in struct_set)
            if write_me:
                for l in curr_block: out.write(l)
            curr_tag = None
            curr_scores = None
            curr_block = None
        elif curr_block is None:
            if line == "\n" or line.startswith("#"): out.write(line)
            else: raise ValueError("Logical failure: non-comment outside of POSE / END_POSE")
        else:
            curr_block.append(line)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
