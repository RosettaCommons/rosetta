#! /usr/bin/env python2
import sys,os

from Types import *
from Molecule import MoleculeClass
from AtomTypeClassifier import AtomTypeClassifier
from TorsionAssigner import TorsionAssigner

MYFILE = os.path.abspath(__file__)
direc = MYFILE.replace('mol2genparams.py','')

def run_mol2(mol2file, torsassigner, report=True, debug=False, prefix=''):
    puckering_as_chi = False
    if '--puckering_chi' in sys.argv: puckering_as_chi = True

    molecule = MoleculeClass(mol2file,verbose=debug)

    classifier = AtomTypeClassifier()
    classifier.apply_to_molecule(molecule)
    classifier.assert_H(molecule)
    
    molecule.assign_rotable_torsions(verbose=debug)

    torsassigner.assign(molecule,debug=debug)

    if report:
        molecule.report_paramsfile('%s.params'%prefix,
                                   report_Hapol_chi=False,
                                   report_puckering_chi=puckering_as_chi,
                                   report_as_atype=True,
                                   )

        molecule.report_pdbfile('%s_0001.pdb'%prefix)

    if '--elec_cp_rep' in sys.argv:
        molecule.report_elec_cp_rep('%s.elec_cp_rep'%prefix)

    if '--elec_grp_def' in sys.argv:
        molecule.report_grpdeffile('%s.grpdef'%prefix)

def main(mol2files):
    torsassigner = TorsionAssigner('%s/torsions.ref'%direc)

    debug = False
    report = True
    prefix = False
    if '--debug' in sys.argv: debug = True
    if '--no_output' in sys.argv: report = False
    if '--prefix' in sys.argv: prefix = sys.argv[sys.argv.index('--prefix')+1]

    for mol2file in mol2files:
        if mol2file.startswith('#'): continue
        print(mol2file)
        try:
            if not prefix:
                prefix_loc = mol2file.split('/')[-1].replace('.mol2','')
                run_mol2(mol2file, torsassigner, report=report, debug=debug, prefix=prefix_loc)
            else:
                run_mol2(mol2file, torsassigner, report=report, debug=debug, prefix=prefix)
        except:
            print('failed on generating params file for ', mol2file)
            if debug:
                run_mol2(mol2file, torsassigner, report=report, debug=True)

if __name__ == "__main__":
    if '-s' not in sys.argv and '-l' not in sys.argv:
        print("usage: python mol2genparams.py [-s mol2file or -l mol2filelist]")
    else:
        if '-s' in sys.argv:
            mol2file = sys.argv[sys.argv.index('-s')+1]
            main([mol2file])
        elif '-l' in sys.argv:
            mol2list = sys.argv[sys.argv.index('-l')+1]
            mol2s = [l[:-1] for l in file(mol2list)]
            main(mol2s)
            
