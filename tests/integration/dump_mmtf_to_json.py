#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
from __future__ import print_function

try:
    import mmtf
except ImportError:
    print("This script requires the Python MMTF parser to be installed. Aborting.")
    exit(0)

import msgpack
import numpy
import sys
import json

def to_json(filename):
    data = {}

    # As of this writing, the mmtf package doesn't extract any of the Properties entries -- do that directly from the msgpack blob.
    with open(filename, 'rb') as fh:
        raw_data = msgpack.load(fh, raw=False)
    for field in ['atomProperties', 'groupProperties', 'chainProperties', 'modelProperties', 'extraProperties', 'bondProperties']:
        if field in raw_data:
            data[field] = raw_data[field]

    struct = mmtf.parse(filename)
    # hopefully all of the relevant data is in the __dict__ element
    # (I'm not sure how to extract it programatically and generally, otherwise.)
    rawdata = struct.__dict__
    # Need to convert certain data to JSON-able format
    for field in rawdata:
        #print( type(rawdata[field]) )
        if type(rawdata[field]) == numpy.ndarray:
            data[ field ] = rawdata[field].tolist()
        else:
            data[ field ] = rawdata[field]

    with open(filename +'.json','w') as f:
        json.dump(data, f, sort_keys=True, indent=2, separators=(',', ': '))

if __name__ == "__main__":
    for filename in sys.argv[1:]:
        to_json(filename)
