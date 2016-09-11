#!/usr/bin/env phenix.python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import os.path
import shutil
import imp

file_path = os.path.split( os.path.abspath(__file__) ) [0]
imp.load_source('erraser_util', file_path + '/erraser_util.py')

from erraser_util import *

#####################################################
print '###################################'
print 'Starting RNA_rosetta_ready_set.py...'
start_time=time.time()
#####Input options#####################################
input_pdb = parse_options(sys.argv, 'pdb', '')
out_pdb = parse_options(sys.argv, 'out_pdb', basename(input_pdb).replace('.pdb', '_ready_set.pdb') )
renumbering = parse_options(sys.argv, 'renumbering', 'False')

if input_pdb =="" :
    error_exit("Error: USER need to specify -pdb option")
check_path_exist(input_pdb)

if exists( out_pdb ) :
    print 'Output pdb %s already exists... Remove it.' % out_pdb
    remove(out_pdb)
#######File Paths#######################################################
base_dir = os.getcwd()
temp_file = '%s/%s_temp_rs.pdb' % (base_dir, basename(input_pdb).replace('.pdb', '') )
input_pdb = abspath(input_pdb)
out_pdb = abspath(out_pdb)

[res_conversion_list, fixed_res_final, cutpoint_final, CRYST1_line] = pdb2rosetta(input_pdb, temp_file, renumbering=renumbering)
if renumbering:
    rna_rosetta_ready_set(temp_file, out_pdb)
else:
    shutil.move(temp_file, out_pdb)
remove(temp_file)

sys.stdout.write( "cutpoint in Rosetta pdb file: " )
for i in cutpoint_final :
    sys.stdout.write( str(i) + ' ' )
sys.stdout.write('\n')

total_time=time.time()-start_time
print '\n', "DONE!...Total time taken= %f seconds" %(total_time) , '\n'
print '###################################'

