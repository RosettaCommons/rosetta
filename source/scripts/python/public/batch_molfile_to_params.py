#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

'''
A wrapper around molfile_to_params.py for generating large numbers of params files and giving them unique names that don't conflict with anything in the standard Rosetta database.  The directory that this script produces can be directly parsed with rosetta using the -in:file:extra_res_batch path

Author: Sam DeLuca
'''
import itertools
import os
import subprocess
import shutil
import fnmatch
import re
import sys
from optparse import OptionParser
from param_utils import *
mol_to_params = "~/rosetta/rosetta_source/src/python/apps/public/molfile_to_params.py"

char_set = ['0','1','2','3','4','5',
            '6','7','8','9','A','B',
            'C','D','E','F','G','H',
            'I','J','K','L','M','N',
            'O','P','Q','R','S','T',
            'U','V','W','X','Y','Z']


def make_params(mol_path,ligand_name,base_name):
    '''Make the params file, append the conformers, and move the output into the proper directory'''
    logfile = "params/"+base_name+"/log.txt"
    log = open(logfile,'a')
    params_cmd = mol_to_params+" -n " +ligand_name + " " + mol_path
    params_child = subprocess.call(params_cmd,shell=True,stdout=log,stderr=log)
    log.close()
    if params_child is not 0:
        return params_child

    pdb_path = "params/"+base_name+"/"+ligand_name+"_conformers.pdb"
    pdb_file = open(pdb_path,'w')
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file,ligand_name+"_*.pdb"):
            conformer =open(file,'r')
            pdb_file.writelines(conformer.readlines())
            conformer.close()
            os.remove(file)

    pdb_file.close()

    if os.path.exists(ligand_name+".params"):
        paramsfile = open(ligand_name+".params",'a')
    else:
        return 1
    paramsfile.write("PDB_ROTAMERS "+ ligand_name+"_conformers.pdb\n")
    paramsfile.close()
    shutil.move(ligand_name+".params", "params/"+base_name+"/"+ligand_name+".params")
    return params_child



if __name__ == "__main__":
    usage = "%prog -d path/to/rosetta/database --script_path=/path/to/molfile_to_params.py molfile_list.txt"
    parser=OptionParser(usage)
    parser.add_option("-d", dest ="database",help="path to minirosetta database",default=1)
    parser.add_option("--script_path",dest="script",help="location of the molfile_to_params script",default = "")
    parser.add_option("--exclusion_list",dest="excluded",help="list of ligand names to manually exclude",default="")
    (options, args) = parser.parse_args()

    if not os.path.exists("params"):
        os.mkdir("params")

    if options.script != "":
        mol_to_params = options.script


    if not os.path.exists(mol_to_params):
        parser.error("ERROR: make sure to specify the path to the molfile_to_params script with the option --params")

    if len(args) != 1:
        parser.error("ERROR: you must specify the path to a file containing a list of molfile paths")

    molfile_list_path = args[0]
    ligand_names = itertools.product(char_set,repeat=3)

    disallowed_ligands = get_disallowed_ligands(options.database+"/")

    if(options.excluded != ""):
        exclude_list = open(options.excluded,'r')
        for line in exclude_list:
            line = line.rstrip()
            disallowed_ligands.add(line)
        exclude_list.close()


    molfile_list = open(molfile_list_path,'r')
    for molfile in molfile_list:
        molfile = molfile.strip()
        while True:
            ligand_name = ligand_names.next()
            ligand_name = "".join(ligand_name)
            if ligand_name not in disallowed_ligands:
                break
        mol_base_name = molfile.split("/").pop().split(".")[0]

        if os.path.exists("params/"+mol_base_name+"/"+ligand_name+".params"):
            print "ligand " + mol_base_name + " already processed, continuing"
            continue
        if not os.path.exists("params/"+mol_base_name):
            os.mkdir("params/"+mol_base_name)
        print ligand_name,mol_base_name
        params_status = make_params(molfile,ligand_name,mol_base_name)
        if params_status is not 0:
            print "WARNING: failed to generate params for " +mol_base_name

