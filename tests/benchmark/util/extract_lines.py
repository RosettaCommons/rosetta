#!/bin/env python

"""extract_lines.py - splits a file based on the the presence/abscence of lines in another file.

Usage

extract_lines.py source.txt lines.txt match.txt nonmatch.txt

For each line in source.tx, if it matches a line in lines.txt, it will be placed into match.txt.
If it doesn't match, it will be placed in nonmatch.txt.
The length of match.txt plus nonmatch.txt should match source.txt.
(That is, `cat match.txt nonmatch.txt | sort` and `sort source.txt` should be identical.)

Note: The current version deletes all numbers from the line before comparison to the reference values
(but outputs with the numbers intact) - this is to remove sensitivity to line number changes.
"""

import sys
import string
import codecs

if len(sys.argv) != 5:
    print __doc__
    exit()

def remove_numbers( s ):
    # return s.translate(None, string.digits) # This doesn't work if s is unicode
    return ''.join([ c for c in s if not c.isdigit() ])

with codecs.open(sys.argv[2], encoding='utf-8', errors='replace') as f:
    ref_lines = set( remove_numbers( line ) for line in f.readlines() )

with codecs.open(sys.argv[1], encoding='utf-8', errors='replace') as f:
    source_lines = f.readlines()

match = []
nonmatch = []

for line in source_lines:
    if remove_numbers(line) in ref_lines:
        match.append( line )
    else:
        nonmatch.append( line )

with codecs.open(sys.argv[3], 'w', encoding='utf-8', errors='replace') as f:
    f.write( ''.join(match) )

with codecs.open(sys.argv[4], 'w', encoding='utf-8', errors='replace') as f:
    f.write( ''.join(nonmatch) )
