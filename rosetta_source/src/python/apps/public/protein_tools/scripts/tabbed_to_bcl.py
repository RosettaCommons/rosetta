#!/usr/bin/env python2.5
from optparse import OptionParser
import sys
from rosettautil.util import fileutil
from rosettautil.bcl import file_formats
usage = "%prog [options] tabbed_input.txt bcl_output.txt"
parser = OptionParser(usage)
(options,args) = parser.parse_args()

in_file = fileutil.universal_open(args[0],"r")
vector_data = file_formats.list_of_2D_vectors()
for line in in_file:
    line = line.split()
    if len(line) == 0:
        continue
    vector_data.add_record(line[0],line[1])
in_file.close()

vector_data.write_bcl_file(args[1])