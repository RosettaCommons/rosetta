#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# Author: Brian D. Weitzner (weitzner@uw.edu)
from __future__ import print_function

from collections import defaultdict
from os import path
import argparse
import gzip
import sys


def _fill_chain_data(fdata, cd):
    {cd[l[21]].append(l) for l in fdata if l.startswith('ATOM')}


def _read_gzipped_pdb_file(fname, cd):
    with gzip.open(fname, 'rb') as f:
        _fill_chain_data([l.decode('utf-8', errors="replace") for l in f.readlines()], cd)


def _read_pdb_file(fname, cd):
    with open(fname, 'r') as f:
        _fill_chain_data(f.readlines(), cd)


def _write_output_file(fname, cd, chains, gzip_results=False):
    open_fxn = gzip.open if gzip_results else open
    with open_fxn(fname, 'w') as f:
        for chain in chains:
            if gzip_results:
                f.writelines([l.encode('utf-8') for l in cd[chain]])
            else:
                f.writelines(cd[chain])


def reorder_chains(in_file, out_file, chains):
    """This function reads a PDB formatted file, reorders the chains as
    requested by the user, and writes the new PDB file to disk.
    """
    fbase, file_extension = path.splitext(in_file)

    open_fxn = _read_pdb_file
    if file_extension == '.gz':
        open_fxn = _read_gzipped_pdb_file
        _, file_extension = path.splitext(fbase)

    if file_extension != '.pdb':
        sys.exit('Input file must be PDB-formatted file. Exiting.')

    chain_data = defaultdict(list)
    open_fxn(in_file, chain_data)

    # we'll iterate over the list of chains twice to separate error checking
    # from file writing while only opening the file one time.
    for chain in chains:
        if chain not in chain_data.keys():
            sys.exit('Chain "{}" not in input. Exiting.'.format(chain))

    _, file_extension = path.splitext(out_file)
    _write_output_file(out_file, chain_data, chains, file_extension == '.gz')


def main(argv):
    """This script is used to reorder the chains in a PDB formatted file.
    NOTE: Only ATOM records are retained.
    """
    parser = argparse.ArgumentParser(description='Program')
    parser.add_argument('-i', '--in_file', action='store', type=str,
                        required=True,
                        help='input PDB file')
    parser.add_argument('-o', '--out_file', action='store', type=str,
                        required=True,
                        help='output PDB file')
    parser.add_argument('-c', '--chains', action='store', type=str,
                        required=True, nargs='+',
                        help='One-letter chainIDs in the desired order')
    args = parser.parse_args()
    reorder_chains(args.in_file, args.out_file, args.chains)


if __name__ == '__main__':
    main(sys.argv)
