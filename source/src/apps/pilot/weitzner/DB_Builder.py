#!/usr/bin/env python
# :noTabs=true: :collapseFolds=1:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   DB_Builder.py
## @brief  Create SQLite DB for RamachandranEnergy2B energey
## @author Sergey Lyskov

import sys
import sqlite3
import optparse

global connection
global connection_cursor

def create_db( db_name ):
    global connection
    global connection_cursor
    
    connection = sqlite3.connect( db_name )
    #connection = sqlite3.connect(':memory:')

    connection_cursor = connection.cursor()

    connection_cursor.execute('''
        CREATE TABLE potential
            (id INTEGER PRIMARY KEY,
            central_aa TEXT, neighbor_aa TEXT,
            left_right INTEGER,
            phi INTEGER, psi INTEGER,
            probability REAL,
            cumulative_probability REAL
            )
    ''')


def add_data(central_aa, neighbor_aa, left_right, phi, psi, probability, cumulative_probability):
    connection_cursor.execute('INSERT INTO potential (id, central_aa, neighbor_aa, left_right, phi, psi, probability, cumulative_probability) VALUES (null, ?, ?, ?, ?, ?, ?, ?)',
              (central_aa, neighbor_aa, bool(left_right), int(phi), (psi), float(probability), float( cumulative_probability ) ) )

def read_data_file( filename ):
    # Open the file, read the lines in and close it.
    try:
        f = file( filename, 'r' )
    except IOError:
        sys.stderr.write( filename + " does not exist. Use the --help flag.\n" )    
        sys.exit()

    lines = f.readlines()
    f.close()
    if len( lines ) == 0:
        sys.stderr.write( filename + " is empty!\n" )   
        sys.exit()
    return lines
    
def main(argv): 
 
    # Missing parameters are handled by providing default values
    usage = "usage: %prog -f input_file -d db_file"
    desc  = "Converts a Dunbrack Neighbor-dependent Ramachandran Distribution file to a sqlite3 database."
    parser = optparse.OptionParser(usage=usage, description=desc)
    
    # Get command line options - use default values if nothing is provided
    parser.add_option("-f","--input_file", default="NDRD_TCBIG.txt",
        dest="input_file",help="Specify the Dunbrack Neighbor-dependent Ramachandran Distribution input file."+\
            "Default is \"NDRD_TCBIG.txt\""
    )
    parser.add_option("-d","--db_file", default="potentials.db",
        dest="db_file",help="Specify the output file. Default is"+\
            " \"potentials.db\""
    )
    # Assign command line flags to variables
    (options, args) = parser.parse_args(args=argv[1:])
    if len(args) != 2:
        sys.stderr.write("Incomplete input given, using default values.\n"+
            "Use \'python DB_Builder.py --help\' to see all options\n")
    input_file = options.input_file
    db_file = options.db_file
    
    create_db( db_file )
    data = read_data_file( input_file )
    
    # Loop through the lines to populate the database
    for datum in data:
        if datum[0] != '#' and len( datum ) > 1:
            datum.rstrip()
            entry = datum.split()
            add_data( entry[ 0 ], entry[ 2 ], int( entry[ 1 ] == 'right' ), entry[ 3 ], entry[ 4 ], entry[ 5 ], entry[ 7 ] )
    connection.commit()
    
if __name__  == "__main__": main(sys.argv)
