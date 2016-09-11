#!/usr/bin/env phenix.python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import os.path
import imp

try :
    import erraser_util 
except :
    file_path = os.path.split( os.path.abspath(__file__) ) [0]
    imp.load_source('erraser_util', file_path + '/erraser_util.py')

from erraser_util import *

#####################################################
print '###################################'
print 'Starting RNA_rosetta_to_pdb.py...'
start_time=time.time()
#####Input options#####################################
input_pdb = parse_options(sys.argv, 'pdb', '')
out_pdb = parse_options(sys.argv, 'out_pdb', basename(input_pdb).replace('.pdb', '_standard.pdb') )

if input_pdb =="" :
    error_exit("Error: USER need to specify -pdb option")
check_path_exist(input_pdb)

if exists( out_pdb ) :
    print 'Output pdb %s already exists... Remove it.' % out_pdb
    remove(out_pdb)

rosetta2std_pdb(input_pdb, 'temp_rs2pdb.pdb')
regularize_OP1_OP2('temp_rs2pdb.pdb', out_pdb)
remove('temp_rs2pdb.pdb')

total_time=time.time()-start_time
print '\n', "DONE!...Total time taken= %f seconds" %(total_time) , '\n'
print '###################################'

