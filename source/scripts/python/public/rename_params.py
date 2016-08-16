#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

'''
Given a directory produced by batch_molfile_to_params.py, rename all the params such that no params
names conflict with anything in the database, on in a list of reserved names

Author: Sam DeLuca
'''

from optparse import OptionParser
from param_utils import *
import glob
import os
import itertools
import shutil

if __name__ == "__main__":

    char_set = ['0','1','2','3','4','5',
                '6','7','8','9','A','B',
                'C','D','E','F','G','H',
                'I','J','K','L','M','N',
                'O','P','Q','R','S','T',
                'U','V','W','X','Y','Z']

    usage = "%prog -d path/to/database params_directory"
    parser=OptionParser(usage)
    parser.add_option("-d", dest ="database",help="path to minirosetta database",default="")
    parser.add_option("--exclusion_list",dest="excluded",help="list of ligand names to manually exclude",default="")
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error("You must specify a params directory")

    if options.database == "":
        parser.error("You must specify a rosetta database with -d")

    param_dir = args[0]

    ligand_names = itertools.product(char_set,repeat=3)
    disallowed_ligands = get_disallowed_ligands(options.database)

    if(options.excluded != ""):
        exclude_list = open(options.excluded,'r')
        for line in exclude_list:
            line = line.rstrip()
            disallowed_ligands.add(line)
        exclude_list.close()

    pdb_path_template = "%(subdir)s/%(param)s.pdb"
    conformer_path_template = "%(subdir)s/%(param)s_conformers.pdb"

    for subdir in glob.glob(param_dir+"/*"):
        if not os.path.isdir(subdir):
            continue

        while True:
            ligand_name = ligand_names.next()
            ligand_name = "".join(ligand_name)
            if ligand_name not in disallowed_ligands:
                break

        params_path = glob.glob(subdir+"/*.params")[0]
        param_name = params_path.split("/")[-1].split(".")[0]
        pdb_path = pdb_path_template % {"subdir" : subdir, "param" : param_name}
        conformer_path = conformer_path_template % {"subdir" : subdir, "param" : param_name}

        new_pdb_path = pdb_path_template % {"subdir" : subdir, "param" : ligand_name}
        new_conformer_path = conformer_path_template % {"subdir" : subdir, "param" : ligand_name}
        new_params_path = subdir+"/"+ligand_name+".params"
        rename_pdb_file(pdb_path,ligand_name)
        rename_pdb_file(conformer_path,ligand_name)
        rename_param_file(params_path,ligand_name,new_conformer_path.split("/")[-1])

        shutil.move(pdb_path,new_pdb_path)
        shutil.move(conformer_path, new_conformer_path)
        shutil.move(params_path,new_params_path)


