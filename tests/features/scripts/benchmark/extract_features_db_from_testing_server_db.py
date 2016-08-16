#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


##################################################################################
# Usage:                                                                         #
#  extract_features_db_from_testing_server_db.py features_db_key output_dir      #
#                                                                                #
#  *feature_database_key: key in features_scientific_benchmark.features_database #
#  *output_dir: where the features_database should be extracted                  #
##################################################################################

import fnmatch, os, re, subprocess, sys

def retrieve_features_db_fname(features_db_key, output_dir):
    stmt_template = '''
SELECT features_db_fname
FROM features_scientific_benchmark.features_database
WHERE key=$(features_db_key);
'''
    stmt = stmt_template % { "features_db_key" : features_db_key }

    ##############################################################
    # Sergey, write code to write stmt to make this work         #
    #                                                            #
    ##############################################################
    features_db_fname_tail = "features_sample_source_id_datecode.db3"

    features_db_fname = os.path.join(output_dir, features_db_fname_tail)
    return features_db_fname

def retrieve_features_db_parts(features_db_key, output_dir):
    stmt_template = '''
SELECT part_fname, part_file_contents
FROM   features_scientific_benchmark.features_database_parts
WHERE  features_db_key=%(features_db_key)s;
'''
    stmt = stmt_template % { "features_db_key" : features_db_key }
    
    ##############################################################
    # Sergey, write code to write stmt to make this work         #
    #                                                            #
    ##############################################################
    part_fnames = []
    for row in []: # query database with select statement
        part_fname = "" # extract part_fname from row
        part_file_contents = None # extract part_file_contents from row
        
        part_fnames.append(part_fname)
        f = open(part_fname, 'wb')
        f.write(part_file_contents)
        f.close()

    return part_fnames


def extract_features_db(features_db_key, output_dir):

    features_db_fname = retrieve_features_db_fname(features_db_key, output_dir)
    part_fnames = retrieve_features_db_parts(features_db_key, output_dir)
    
    p = subprocess.Popen(
        args = ["cat"] + part_fnames + [">", features_db_fname+".tar.gz"],
        shell = True)
    p.wait()
    [ os.remove(part_fnames) for part_fname in part_fnames ]
    
    p = subprocess.Popen(["tar", "-xzf", features_db_fname+".tar.gz"])
    p.wait()
    os.remove(features_db_fname+".tar.gz")

    return features_db_fname
    

if __name__ == "__main__":
    features_db_key = sys.argv[1]
    output_dir = sys.argv[2] 
    extract_features_db(features_db_key, output_dir)
