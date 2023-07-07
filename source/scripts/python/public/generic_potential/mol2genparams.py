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
import tarfile
from ReportTorsions import report_raw_torsions
from utils import GetMol2Block, validate_restype
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
        prefix = os.path.basename(mol2file).replace('.mol2','')
    else:
        prefix = option.opt.prefix

    molecule = MoleculeClass(mol2file,option)
    if option.opt.outdir:
        outdir = option.opt.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    else:
        outdir="./"

    if not option.opt.no_output:
        params_outfn = os.path.join(outdir, '%s.params'%prefix)
        molecule.report_paramsfile(params_outfn)
        if not option.opt.no_pdb:
            pdb_outfn = os.path.join(outdir, '%s_0001.pdb'%prefix)
            molecule.report_pdbfile(pdb_outfn)

    if option.opt.report_funcgrp:
        classifier = FunctionalGroupClassifier()
        classifier.apply_to_molecule(molecule)
        molecule.report_functional_grps(sys.stdout)
        
    if option.opt.write_elec_cp_rep:
        molecule.report_elec_cp_rep('%s.elec_cp_rep'%prefix)

    if option.opt.write_elec_grpdef:
        molecule.report_grpdeffile('%s.grpdef'%prefix)
        
    if option.opt.write_raw_torsion:
        torsion_outfn = os.path.join(outdir, '%s.torsion'%prefix)
        report_raw_torsions( molecule, torsion_outfn )

def run_mol2_block(mol2file, option):
    with open(mol2file, 'r') as infh:
        for mol2block in GetMol2Block(infh):
            try:
                molecule = MoleculeClass(mol2block, option)
            except:
                print("Failed:", mol2block.split('\n')[1])
                if option.opt.debug:
                    print(mol2block)
                    molecule = MoleculeClass(mol2block, option)
                    input("Press Enter to continue...")
                continue
            prefix = molecule.name
            if not validate_restype(prefix):
                print(f"{prefix} is not a valid restype name")
                continue
            if option.opt.outdir:
                os.makedirs(option.opt.outdir, exist_ok=True)
                outdir = option.opt.outdir
            else:
                outdir="./"

            if not option.opt.no_output:
                params_outfn = os.path.join(outdir, '%s.params'%prefix)
                molecule.report_paramsfile(params_outfn)
                if not option.opt.no_pdb:
                    pdb_outfn = os.path.join(outdir, '%s_0001.pdb'%prefix)
                    molecule.report_pdbfile(pdb_outfn)

            if option.opt.report_funcgrp:
                classifier = FunctionalGroupClassifier()
                classifier.apply_to_molecule(molecule)
                molecule.report_functional_grps(sys.stdout)
                
            if option.opt.write_elec_cp_rep:
                molecule.report_elec_cp_rep('%s.elec_cp_rep'%prefix)

            if option.opt.write_elec_grpdef:
                molecule.report_grpdeffile('%s.grpdef'%prefix)
                
            if option.opt.write_raw_torsion:
                torsion_outfn = os.path.join(outdir, '%s.torsion'%prefix)
                report_raw_torsions( molecule, torsion_outfn )


def run_mol2_tar_gz(mol2file,option):
    mode="r"
        
    with tarfile.open(mol2file, mode) as intar:
        c = 0
        for member in intar.getmembers():
            name = member.name
            print(name)
            if os.path.basename(name).startswith("._"): continue
            if not name.endswith(".mol2"):
                print("Warning: skipping %s due to non-mol2 extension."%os.path.basename(name))
                continue
            
            fobj = intar.extractfile(member)
            if fobj is None:
                continue
            c += 1
            if option.opt.prefix == None:
                prefix = os.path.basename(name).split(".")[0]
            else:
                prefix = option.opt.prefix
            molecule = MoleculeClass(name,option,fobj)
            if not molecule.stat:
                print(f"{name} failed reading.")
                continue
            
            if option.opt.outdir:
                outdir = option.opt.outdir
                os.makedirs(outdir, exist_ok=True)
            else:
                outdir = "./"

            if not option.opt.no_output:
                params_outfn = os.path.join(outdir, '%s.params'%prefix)
                molecule.report_paramsfile(params_outfn)
                if not option.opt.no_pdb:
                    pdb_outfn = os.path.join(outdir, '%s_0001.pdb'%prefix)
                    molecule.report_pdbfile(pdb_outfn)

            if option.opt.report_funcgrp:
                classifier = FunctionalGroupClassifier()
                classifier.apply_to_molecule(molecule)
                molecule.report_functional_grps(sys.stdout)
        
            if option.opt.write_elec_cp_rep:
                molecule.report_elec_cp_rep('%s.elec_cp_rep'%prefix)

            if option.opt.write_elec_grpdef:
                molecule.report_grpdeffile('%s.grpdef'%prefix)
            
            if option.opt.write_raw_torsion:
                torsion_outfn = os.path.join(outdir, 'torsions', '%s.torsion'%prefix)
                report_raw_torsions( molecule, torsion_outfn )

    print(f"Attempted {c} entries in {mol2file}")

def main(option):
    mol2files = option.opt.inputs
    
    for mol2file in mol2files:
        if mol2file.startswith('#'): continue
        print(mol2file)
        if not os.path.exists(mol2file):
            print("WARNING: Cannot find mol2file %s, skip!"%mol2file)
            continue
        try:
            if mol2file.endswith(".mol2"):
                if option.opt.multimol2:
                    run_mol2_block(mol2file, option)
                else:
                    run_mol2(mol2file,option)
            elif mol2file.endswith(".gz") or mol2file.endswith(".tar"):
                run_mol2_tar_gz(mol2file,option)
        except:
            print('failed on generating params file for ', mol2file)
            if option.opt.debug:
                if mol2file.endswith(".mol2"):
                    run_mol2(mol2file,option)
                elif mol2file.endswith(".gz") or mol2file.endswith(".tar"):
                    run_mol2_tar_gz(mol2file,option)
                
        option.resname_counter += 1 #add counter

if __name__ == "__main__":
    option = OptionClass(sys.argv)

    main(option)
            
