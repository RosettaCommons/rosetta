import sys,os

from Types import *
from Molecule import MoleculeClass
from AtomTypeClassifier import AtomTypeClassifier
from TorsionAssigner import TorsionAssigner
from utils import *


def report_raw_torsions(molecule, outfile):
    # BOND_ORDERS = [1, # single
    #        2, # double
    #        3, # triple
    #        4, # aromatic
    #        5, # orbital, not used
    #        9, # ring, not used; will be used once Rosetta machinery gets updated
    #        ]
    content = ["#torsion type, bond order, torsion angle\n"]
    for torsion in molecule.torsions:
        [i0,i1,i2,i3] = torsion
    
        at0 = ACLASS_ID[molecule.atms[i0].aclass]
        at1 = ACLASS_ID[molecule.atms[i1].aclass]
        at2 = ACLASS_ID[molecule.atms[i2].aclass]
        at3 = ACLASS_ID[molecule.atms[i3].aclass]
    
        xyz0 = molecule.xyz[i0]
        xyz1 = molecule.xyz[i1]
        xyz2 = molecule.xyz[i2]
        xyz3 = molecule.xyz[i3]
    
        t = dihedral(xyz0,xyz1,xyz2,xyz3)
        
        bt = molecule.bond_order(i1,i2)
        #print molecule.atms[i1],molecule.atms[i2]
    
        content.append("%6s %6s %6s %6s %2d %.2f \n"%(at0, at1, at2, at3, bt, t))
    dirname = os.path.dirname(outfile)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    with open(outfile,'w') as output:
        output.writelines(content)
    print("Saved: %s"%outfile)
