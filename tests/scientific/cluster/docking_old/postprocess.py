#!/usr/bin/env python
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# :noTabs=true:

import sys


#
# Convert string to float if possible, convert 'NA' string to Python None value
#
def convert(v):
    if v == 'NA': return None

    try: v = float(v)
    except (ValueError):
        pass

    return v



def main(argv):
    '''
    Extract docking test result and convert them to YAML, saving to '.results.yaml'.
    '''
    data = file(argv[0]).read()

    lines = data.split('\n')
    legend = lines[0].split(';')
    #print 'lines', lines
    #print 'legend', legend

    res = {}
    for l in lines[1:]:
        d = l.split(';')
        if len(d) < len(legend): continue
        filename = d[0]
        res[filename] = {}
        for i in range(1, len(legend)):
            #print i, filename,  d[i]
            value = convert( d[i] )
            res[ filename ][ legend[i] ] = value
    #print res

    file('.results.yaml', 'w').write( str(res) )





if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
