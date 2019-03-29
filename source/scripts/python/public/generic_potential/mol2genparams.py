#!/usr/bin/env python3
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
'''
Functions and executable for taking a ligand from an MDL Molfile
and writing .params files with generalized atom types in Rosetta.
See main() for usage or run with --help.

Author: Hahnbeom Park and Frank DiMaio 
'''
from __future__ import print_function

import sys,os

from Types import *
from BasicClasses import OptionClass
from Molecule import MoleculeClass
from AtomTypeClassifier import FunctionalGroupClassifier

#### check if libraries are installed
import importlib

die_if_not_exist = ['numpy','scipy']
for lib in die_if_not_exist:
    die = (importlib.util.find_spec(lib) is None)
    if die:
        sys.exit('this script requires %s! die.'%lib)
    
############

MYFILE = os.path.abspath(__file__)
direc = MYFILE.replace('mol2genparams.py','')

def run_mol2(mol2file,option):
    # local prefix for the outputs
    if option.opt.prefix == None:
        prefix = mol2file.split('/')[-1].replace('.mol2','')
    else:
        prefix = option.opt.prefix

    molecule = MoleculeClass(mol2file,option)

    if not option.opt.no_output:
        molecule.report_paramsfile('%s.params'%prefix)
        molecule.report_pdbfile('%s_0001.pdb'%prefix)

    if option.opt.report_funcgrp:
        classifier = FunctionalGroupClassifier()
        classifier.apply_to_molecule(molecule)
        molecule.report_functional_grps(sys.stdout)
        
    if option.opt.write_elec_cp_rep:
        molecule.report_elec_cp_rep('%s.elec_cp_rep'%prefix)

    if option.opt.write_elec_grpdef:
        molecule.report_grpdeffile('%s.grpdef'%prefix)

def main(option):
    mol2files = option.opt.inputs
    
    for mol2file in mol2files:
        if mol2file.startswith('#'): continue
        print(mol2file)
        if not os.path.exists(mol2file):
            print("WARNING: Cannot find mol2file %s, skip!"%mol2file)
            continue
        try:
            run_mol2(mol2file,option)
        except:
            print('failed on generating params file for ', mol2file)
            if option.opt.debug:
                run_mol2(mol2file,option)
                
        option.resname_counter += 1 #add counter

if __name__ == "__main__":
    option = OptionClass(sys.argv)

    main(option)
            
