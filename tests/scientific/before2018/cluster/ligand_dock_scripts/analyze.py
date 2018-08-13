#!/usr/bin/env python
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

from sys import argv, exit
from glob import glob
from optparse import OptionParser
from targets import targets
import yaml

def parse_lines(inp):
    cols = None  # mapping of score names to column numbers
    decoys = []  # list of lists of scores values, indexed by cols
    for line in inp:
        if not line.startswith("SCORES "): continue
        if "is_reference_pose 1" in line: continue # skip ref. structures
        f = line.rstrip().split()
        if cols is None:
            colnames = f[::2] # entries 0, 2, 4, ...
            cols = dict(zip(colnames, range(len(colnames))))
        values = f[1::2] # entries 1, 3, 5, ...
        decoys.append([ float(values[cols["total_score"]]), float(values[cols["interface_delta_X"]]), float(values[cols["ligand_rms_no_super_X"]]) ])
    n_total_decoys = len(decoys)
    decoys.sort(key=lambda x: x[0]) # sort by total score
    decoys = decoys[:int(round(0.05*len(decoys)))] # retain best 5%
    decoys.sort(key=lambda x: x[1]) # re-sort by interface delta
    rms_of_first = decoys[0][2]
    rank_below_2 = 0 # zero-based rank, needs +1 for "normal" interpretation
    rms_below_2 = None
    while rank_below_2 < len(decoys):
        if decoys[rank_below_2][2] <= 2.0:
            rms_below_2 = decoys[rank_below_2][2]
            break
        rank_below_2 += 1
    else:
        rank_below_2 = None

    return (rms_of_first, rank_below_2, rms_below_2, n_total_decoys)

def main(argv):
    '''
    Extract SCORE information from ligand atom_tree_diff output files
    and write plain-text and YAML reports for scientific benchmarks.
    '''
    parser = OptionParser(usage="usage: %prog atom_tree_diff1.out atom_tree_diff2.out ...")
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
    sta_result = open("staResult", "w")
    sta_result.write("%-24s%16s%16s%16s\n" % ("FILENAME", "RMS OF #1", "RANK BELOW 2A", "RMS BELOW 2A"))
    
    records = {}
    records['NumOfTargetsWithRMSBelow2'] = 0 

    for target in targets:
      files= glob('output'+'/'+target+'*'+'atom_tree_diff'+'*')
      if len(files) == 0: continue
      lines=[]
      for file in files:
        lines += open(file).readlines()
      (rms_of_first, rank_below_2, rms_below_2, n_total_decoys) = parse_lines(lines)
      if rank_below_2 is not None: rank_below_2 += 1
      #print target, rms_of_first, rank_below_2, rms_below_2
      # For plain text output:
      if rank_below_2 is not None:
          sta_result.write("%-24s%16.2f%16i%16.2f\n" % (target, rms_of_first, rank_below_2, rms_below_2))
      else:
          sta_result.write("%-24s%16.2f%16s%16s\n" % (target, rms_of_first, "-", "-"))
      # For YAML output:
      records[target] = {
          "rms_of_first": rms_of_first,
          "rank_below_2": rank_below_2,
          "rms_below_2": rms_below_2,
      }
      if rms_of_first < 2: records['NumOfTargetsWithRMSBelow2'] += 1
      
    yaml.dump(records, open(".results.yaml", "w"))

    return 0

if __name__ == "__main__":
    exit(main(argv[1:]))
