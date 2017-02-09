#!/usr/bin/env python
##
## (c) Copyright Rosetta Commons Member Institutions.
## (c) This file is part of the Rosetta software suite and is made available under license.
## (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
## (c) For more information, see http://www.rosettacommons.org. Questions about this can be
## (c) addressed to University of Washington UW CoMotion, email: license@u.washington.edu.
##
## @author Rocco Moretti (rmorettiase@gmail.com)

''' select_columns.py - Select only certain columns from a file, given their contents.

Particularly useful for Rosetta scorefiles, as it allows you to select columns by headings.
'''

import sys, os
import argparse
import itertools
import re

def main(args):

    lines = []
    if args.start is not None:
        useline = False
    else:
        useline = True

    with open(args.infile) as f:
        for line in f:
            if not useline:
                if args.start is not None and args.start in line:
                    useline=True
                continue
            if args.end is not None and args.end in line:
                break
            lines.append( line.split(args.delimiter) )

    #transpose
    columns = map(list, itertools.izip_longest(*lines, fillvalue='') )

    selected = []

    #TODO: Fix column selection logic.
    for entry in args.columns:

        if args.numeric:
            colnum = None
            try:
                colnum = int(entry)
            except ValueError:
                colnum = None
            if colnum is not None:
                if colnum <= 0:
                    if args.verbose:
                        sys.stderr.write("Zero-based index selection for "+str(colnum)+'\n')
                    # Allow negative "from the end" indexing
                    selected.append( columns[colnum] )
                    continue
                else:
                    if args.verbose:
                        sys.stderr.write("One-based index selection for "+str(colnum)+'\n')
                    # Convert from 1-based indexing
                    selected.append( columns[colnum-1] )
                    continue
            #Fall-back to text based selection is deliberate

        if args.regex:
            rx = re.compile(entry)

        for column in columns:
            if args.regex:
                if any( rx.search( i ) for i in column ):
                    if args.verbose:
                        sys.stderr.write("Regex-based selection for "+entry+'\n')
                    selected.append( column )
            else:
                # looking for literals
                if entry in column:
                    if args.verbose:
                        sys.stderr.write("Text selection for "+entry+'\n')
                    selected.append(column)

    #OUTPUT

    delim = args.delimiter
    if delim is None:
        #Use spaces to align columns
        for colnum, col in enumerate(selected):
            maxlen = max( len(c) for c in col )
            selected[colnum] = [c.rjust(maxlen) for c in col]
        delim = '  '

    #Back transpose
    lines = map(list, itertools.izip_longest(*selected, fillvalue='') )

    for line in lines:
        if not any(l.strip() for l in line):
            # Ignore completely blank lines
            continue
        print delim.join(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
        epilog='''Usage note: Keep in mind that leading/trailing whitespace is significant with non-default delimiters. You may wish to combine `-d` with `-r`.''')

    parser.add_argument("infile", help="The tabular input file from which you're selecting the columns.")
    parser.add_argument("columns", nargs='+', help="The columns which you want to select.")
    parser.add_argument('-d','--delimiter', default=None, help="Delimeter to use to split columns (default whitespace)")
    parser.add_argument('-n','--numeric', action='store_true', help="Select columns by numbers (for ones which are integers) rather than strings.")
    parser.add_argument('-r','--regex', action='store_true', help="Column designations are interpreted as Python regexes, rather than a literals.")
    parser.add_argument('--start', default=None, help="Don't start until the line *after* the line containing the given text")
    parser.add_argument('--end', default=None, help="End on the line *before* the line containing the given text")
    parser.add_argument('-v','--verbose', action='store_true', help="Print progress messages to stderr, to help with debugging.")

    args = parser.parse_args()
    main(args)

