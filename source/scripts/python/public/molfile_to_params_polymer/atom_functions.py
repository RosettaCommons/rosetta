#!/usr/bin/env python
from __future__ import print_function

import os, sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

# Magic spell to make sure Rosetta python libs are on the PYTHONPATH:
#sys.path.append( os.path.abspath( sys.path[0] ) )
#sys.path.append(os.path.dirname(os.path.abspath(__file__)))

#from rosetta_py.io.mdl_molfile import *
#from rosetta_py.utility.rankorder import argmin
#from rosetta_py.utility import r3

from molfile_to_params_polymer.bond_functions import *
def add_fields_to_atoms(atoms):
    '''Adds a bunch of member variable placeholders that we use.'''
    for atom in atoms:
        atom.orig_name = atom.name # for kinemage output
        atom.pdb_name = ""         # PDB style atom name
        atom.ros_type = ""         # Rosetta atom type
        atom.mm_type = ""          # Molec. mechan. atom type
        atom.is_virtual = False    # for atom typing and charge assignment
        atom.rigid_id = 0          # non-zero id for regions with no rotatable bonds; may span fragments
        atom.fragment_id = 0       # non-zero id for fragments after bond breaking
        atom.conn_bonds = []       # list of cross-fragment bonds to this atom
        atom.is_root = False       # for atom tree
        atom.parent = None         # for atom tree
        atom.children = []         # for atom tree
        atom.stub1 = None          # for internal coords, derived from atom tree
        atom.stub2 = None          # for internal coords, derived from atom tree
        atom.stub3 = None          # for internal coords, derived from atom tree
        atom.input_stub1 = None    # for internal coords, derived from atom tree
        atom.input_stub2 = None    # for internal coords, derived from atom tree
        atom.input_stub3 = None    # for internal coords, derived from atom tree
        atom.d = 0.0               # distance to input_stub1
        atom.theta = 0.0           # 180 - angle with input_stub2 (degrees)
        atom.phi = 0.0             # dihedral from input_stub3 (degrees)
        atom.poly_upper = False    # is upper connect atom for polymer residue type
        atom.poly_lower = False    # is lower connect atom for polymer residue type
        atom.poly_n_bb = False     # is backbone nitrogen for polymer residue type
        atom.poly_ca_bb = False    # is backbone alpha carbon for polymer residue type
        atom.poly_c_bb = False     # is backbone carbonyl carbon for polymer residue type
        atom.poly_o_bb = False     # is backbone carbonyl oxygen for polymer residue type
        atom.poly_backbone = False # convience boolean
        atom.poly_ignore = False   # convience boolean
        atom.pdb_prefix_num = " "  # blah
        atom.pdb_elem = atom.elem  # blah
        atom.pdb_greek_dist = " "  # blah
        atom.pdb_postfix_num = " " # blah

def find_virtual_atoms(atoms):
    '''Atoms whose names start with "V" are virtual, used in enzyme design, etc.'''
    for atom in atoms:
        if atom.name.startswith("V") or atom.name.startswith("X"):
            atom.is_virtual = True
            atom.elem = "X" # "V" is a real element (Vanadium)

def uniquify_atom_names(atoms):
    '''If atom names are not unique, rename/number them.'''
    # First try to pad names to match PDB convention:
    for atom in atoms:
        # If two chars, assume it's an element name
        #if len(atom.name) >= 2 and atom.name[1].isalpha(): atom.name = "%-4s" % atom.name
        #else: atom.name = " %-3s" % atom.name
        # Assume that we got the element name right, and that the atom name starts with the element symbol
        if len(atom.elem) == 1 and len(atom.name) <= 3: atom.name = " %-3s" % atom.name
        else: atom.name = "%-4s" % atom.name
    duplicate_names = False
    atom_names = set()
    for atom in atoms:
        if atom.name in atom_names:
            duplicate_names = True
            break
        atom_names.add(atom.name)
    if not duplicate_names: return
    print( "Atom names contain duplications -- renaming all atoms." )
    # This is potentially bad (> 4 char names) for > 100 atoms:
    #for i,atom in enumerate(atoms):
    #    atom.name = "%s%i" % (atom.elem, i+1)
    # So instead, number each element separately:
    atom_names.clear()
    for atom in atoms:
        i = 1
        while True:
            name = "%2s%i" % (atom.elem, i)
            if name not in atom_names: break
            i += 1
        atom.name = name
        atom_names.add(name)

def assign_rosetta_types(atoms):
    '''Assigns Rosetta atom types.
    Based on Rosetta++ ligand_ns.cc set_rosetta_atom_types().
    This has been tested against the assignments produced by (Jens? Molecule.exe?)
    for the Meiler and Baker 2006 cross docking test set;
    cases where they disagree are usually due to weird
    bond orders (probably mistakes) in the .mol files.

    As I look through these rules, I see that they contain an ad-hoc attempt
    at aromaticity perception.  Not sure how general this is, though.
    '''
    # Helper function, count bonded atoms that match predicate
    # Predicate takes one arg, an atom
    def count_bonded(atom, pred):
        return sum(1 for bond in atom.bonds if not bond.a2.is_virtual and pred(bond.a2))
    # A predicate for count_bonded()
    def is_aromatic(atom):
        return any( bond.order == Bond.AROMATIC for bond in atom.bonds if not bond.a2.is_virtual )
    # For carbon and nitrogen, is_saturated implies SP3 hybridization
    # Not quite true:  look at ring N in 1aq1 -- 3 "single" bonds but flat (like Trp)
    def is_saturated(atom):
        return all( bond.order == Bond.SINGLE for bond in atom.bonds if not bond.a2.is_virtual )
    # For each atom, analyze bonding pattern to determine type
    for i, a in enumerate(atoms): # i just used for debugging output
        # H, C, O, N have complicated rules.
        # Everything else maps to a single atom type.
        if a.is_virtual:
            a.ros_type = "VIRT"
        elif a.poly_n_bb:
            heavy_nbrs = count_bonded(a, lambda x: not x.is_H)
            if heavy_nbrs >= 3:
                a.ros_type = "Npro"
            else:
                a.ros_type = "Nbb "
        elif a.poly_ca_bb:  a.ros_type = "CAbb"
        elif a.poly_o_bb:   a.ros_type = "OCbb"
        elif a.poly_c_bb:   a.ros_type = "CObb"
        elif "H"  == a.elem:
            num_aro_C = count_bonded(a, lambda x: (x.elem == "C" and is_aromatic(x)) or x.ros_type == "aroC")
            num_NO = count_bonded(a, lambda x: x.elem == "N" or x.elem == "O")
            num_S = count_bonded(a, lambda x: x.elem == "S")
            
            if any(bond.a2.poly_n_bb for bond in a.bonds): a.ros_type = "HNbb"
            elif num_S >= 1:    a.ros_type = "HS  "
            elif num_NO >= 1:   a.ros_type = "Hpol"
            elif num_aro_C >=1: a.ros_type = "Haro"
            else:               a.ros_type = "Hapo"
        elif "C"  == a.elem:
            num_H = count_bonded(a, lambda x: x.is_H)
            if is_saturated(a):
                if num_H >= 3:      a.ros_type = "CH3 "
                elif num_H == 2:    a.ros_type = "CH2 "
                else:               a.ros_type = "CH1 "
            else:
                num_dbl_nonO = 0; num_aro_nonO = 0; num_aro_N = 0; num_sgl_N = 0;
                a_bonds = [b for b in a.bonds if not b.a2.is_virtual]
                for bond in a_bonds:
                    if bond.order == Bond.SINGLE:
                        if bond.a2.elem == "N": num_sgl_N += 1
                    elif bond.order == Bond.DOUBLE:
                        if bond.a2.elem != "O": num_dbl_nonO += 1
                    elif bond.order == Bond.AROMATIC:
                        if bond.a2.elem != "O": num_aro_nonO += 1
                        if bond.a2.elem == "N": num_aro_N += 1 # really if, not elif
                #print( i+1, a.name, num_aro_nonO, num_dbl_nonO, num_aro_N )
                if num_aro_nonO >= 2: 
                    if num_H > 0:       a.ros_type = "aroC"
                    else:               a.ros_type = "CH0 "
                elif num_dbl_nonO >= 1: 
                    if num_H > 0:       a.ros_type = "aroC"
                    else:               a.ros_type = "CH0 "
                elif num_aro_N >= 1:    a.ros_type = "CNH2"
                elif num_sgl_N == 1:    a.ros_type = "CNH2"
                else:                   a.ros_type = "COO "
        elif "N"  == a.elem:
            num_H = count_bonded(a, lambda x: x.is_H)
            heavy_nbrs = count_bonded(a, lambda x: not x.is_H)
            assert num_H + heavy_nbrs == len(a.bonds)
            if num_H >= 3:              a.ros_type = "Nlys" # carries a VERY high desolvation penalty
            # Not totally sure about this one, may want Ntrp instead if more than one heavy neighbor:
            elif num_H == 2:            a.ros_type = "NH2O" # Narg would also be a possibility, but they're fairly similar
            elif num_H == 1:
                if heavy_nbrs <= 2:     a.ros_type = "Ntrp" # should always be 2 neighbors, not less
                else:                   a.ros_type = "Ntrp" # Npro? protonated tertiary amine
            else: # num_H == 0
                if heavy_nbrs <= 2:     a.ros_type = "Nhis"
                elif heavy_nbrs == 3:
                    if is_saturated(a): a.ros_type = "Nhis" # deprotonated tertiary amine; need an sp3 hybrid H-bond acceptor type...
                    # This also catches nitro groups -- is that what we want here?
                    else:               a.ros_type = "Npro" # X=[N+](X)X, including nitro groups
                else:                   a.ros_type = "Npro" # quaternary amine
        elif "O"  == a.elem:
            num_H = count_bonded(a, lambda x: x.is_H)
            num_bonds = count_bonded(a, lambda x: True) #len(a.bonds)  using count_bonded() avoid counting virtual bonds/atoms
            bonded_to_N = count_bonded(a, lambda x: x.elem == "N")
            bonded_to_C_to_N = count_bonded(a, lambda x: x.elem == "C" and count_bonded(x, lambda y: y.elem == "N") > 0)
            if is_saturated(a):
                if num_bonds >= 2:
                    unsat_nbrs = count_bonded(a, lambda x: not is_saturated(x))
                    if num_H > 0:   a.ros_type = "OH  " # catches C(=O)OH (Kekule form)
                    elif a.is_ring and a.ring_size < 5:
                                    a.ros_type = "OH  " # small, strained rings leave the O more exposed? (IWD, see 1p8d)
                    elif a.is_ring and unsat_nbrs > 0:
                                    a.ros_type = "Oaro" # catches aromatic O in furan-like rings, though I rarely see these H-bond (IWD)
                    else:           a.ros_type = "OH  " # catches ethers, ROR (IWD, see comment)
                    # The lone pairs on ethers are capable of H-bonding in the same way that alcohols are.
                    # While alkyl ethers are quite non-polar, many others seem to make Hbonds,
                    # such as those attached to phosphates (R-O-PO3), methyls (R-O-CH3), and aromatic rings (R-O-Ph).
                    # It is unclear from the literature how strong these are, and is probably very situation dependent.
                else:               a.ros_type = "OOC " # catches C(=O)[O-] (Kekule form) -- new rule by IWD
            elif num_H > 0:         a.ros_type = "OH  " # catches c(o)oH (aromatic bonds to both O)
            elif bonded_to_N:       a.ros_type = "ONH2"
            # This is a non-standard rule introduced by IWD, agreed to by KWK:
            elif bonded_to_C_to_N:  a.ros_type = "ONH2"
            else:                   a.ros_type = "OOC "
        elif "S"  == a.elem: 
            num_H = count_bonded(a, lambda x: x.is_H)
            if num_H == 1:   a.ros_type = "SH1 "
            else:            a.ros_type = "S   "
        elif "P"  == a.elem: a.ros_type = "Phos"
        elif "F"  == a.elem: a.ros_type = "F   "
        elif "CL" == a.elem: a.ros_type = "Cl  "
        elif "BR" == a.elem: a.ros_type = "Br  "
        elif "I"  == a.elem: a.ros_type = "I   "
        elif "NA" == a.elem: a.ros_type = "Na1p"
        elif "K"  == a.elem: a.ros_type = "K1p "
        elif "MG" == a.elem: a.ros_type = "Mg2p"
        elif "FE" == a.elem: a.ros_type = "Fe3p"
        elif "CA" == a.elem: a.ros_type = "Ca2p"
        elif "ZN" == a.elem: a.ros_type = "Zn2p"
        elif "B" == a.elem: a.ros_type = "Bsp2"
        elif "SE" == a.elem: a.ros_type = "Se  "
        else: raise ValueError("Unknown element '%s'" % a.elem)

def assign_mm_types(atoms, peptoid):
    '''Written by Doug Renfrew strongly influenced by the above function. This _TRYS_ to fill in the CHARMM27 atom types.
    It is rather conservative. It may but X in places in which case you will be tasked with finding the correct type, sorry.
    Rather than the crazy if/else tree (probably more appropriately forest) above I tried to split each type up in to
    its own function. There is unfortunatly some order dependence in calling all the functions but I have tried to avoid it. '''
    # helper functions
    def count_bonded(atom, pred):
        return sum(1 for bond in atom.bonds if not bond.a2.is_virtual and pred(bond.a2))
    def is_aromatic(atom):
        return any( bond.order == Bond.AROMATIC for bond in atom.bonds if not bond.a2.is_virtual )
    def is_saturated(atom):
        return all( bond.order == Bond.SINGLE for bond in atom.bonds if not bond.a2.is_virtual )
    # ----- HYDROGEN -----
    # ----- HYDROGEN -----
    # ----- HYDROGEN -----
    def is_charmm_H(atom, attached):
        ''' For polar hydogen: attached to O or one of 2 H attached to a N '''
        if attached.elem == "O": return True
        elif attached.elem == "N":
            attached_num_H =  count_bonded(attached, lambda x: x.is_H)
            if attached_num_H <= 2: return True
        else: return False
    def is_charmm_HC(atom, attached):
        ''' For nterm hydogen: one of 3 H attached to a N '''
        if attached.elem == "N":
            attached_num_H = count_bonded(attached, lambda x: x.is_H)
            if attached_num_H == 3: return True
            elif is_aromatic(attached) and attached_num_H == 2: return True
        else: return False
    def is_charmm_HA(atom, attached):
        ''' For nonpolar hydrogen: an H attached to a saturated carbon '''
        attached_sat = count_bonded(atom, lambda x: (x.elem == "C" and is_saturated(x)))
        if attached_sat: return True
        else: return False
    def is_charmm_HP(atom, attached):
        ''' For aromatic hydrogen: hydrogen attached to an aromatic carbon '''
        num_aro_C = count_bonded(atom, lambda x: (x.elem == "C" and is_aromatic(x)))
        if num_aro_C >= 1: return True
    def is_charmm_HB(atom, attached, peptoid):
        ''' For backbone hydrogen: hydogen attached to a backbone carbon alpha '''
        if attached.poly_ca_bb and not peptoid: return True
        else: return False
    def is_charmm_HR1(atom, attached):
        ''' For nutral HIS HE1 hydrogen: (A) attached to an aromatic carbon that is between two N and only
        one is prot or (B) attached to and aromatic carbon between an N and C and both ring N are protonated '''
        if is_aromatic(attached) and attached.elem == "C":
            attached_N = [bond.a2 for bond in attached.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "N"]
            attached_C = [bond.a2 for bond in attached.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
            # case A
            if len(attached_N) == 2:
                num_prot_N = sum(1 for an in attached_N if len(an.bonds) == 3 ) # catches H and CH3
                if num_prot_N == 1: return True
                else: return False
            # case B
            elif len(attached_N) == 1 and len(attached_C) == 1:
                attached_C_N = [bond.a2 for bond in attached_C[0].bonds if bond.order == Bond.AROMATIC and  bond.a2.elem == "N"]
                num_sat_N =  sum(1 for an  in attached_N   if len(an.bonds ) == 3 ) # catches H and CH3
                num_sat_N += sum(1 for acn in attached_C_N if len(acn.bonds) == 3 ) # catches H and CH3
                if num_sat_N == 2: return True
                else: return False
            else: return False
        else:return False
    def is_charmm_HR2(atom, attached):
        ''' For protonated HIS HE1: attached to aromatic C b/w two protonated N '''
        if is_aromatic(attached) and attached.elem == "C":
            attached_N = [bond.a2 for bond in attached.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "N"]
            if len(attached_N) == 2:
                num_prot_N = sum(1 for an in attached_N if len(an.bonds) == 3 ) # catches H and CH3
                if num_prot_N == 2: return True
                else: return False
            else: return False
        else:return False
    def is_charmm_HR3(atom, attached):
        ''' For nutral HIS HD2 hydrogen: attached to aromatic C b/w an N and C and one N on ring is protonated'''
        if is_aromatic(attached) and attached.elem == "C":
            attached_N = [bond.a2 for bond in attached.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "N"]
            attached_C = [bond.a2 for bond in attached.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
            if len(attached_N) == 2: return False
            elif len(attached_N) == 1 and len(attached_C) == 1:
                attached_C_N = [bond.a2 for bond in attached_C[0].bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "N"]
                if len(attached_C_N) == 0: return False #The residue is not HIS
                num_sat_N =  sum(1 for an  in attached_N   if len(an.bonds ) == 3 ) # catches H and CH3
                num_sat_N += sum(1 for acn in attached_C_N if len(acn.bonds) == 3 ) # catches H and CH3
                if num_sat_N == 1: return True
                else: return False
            else: return False
        else:return False
    def is_charmm_HS(atom, attached):
        ''' For thiol hydrogen: hydrogen attached to a sulfer '''
        if attached.elem == "S": return True
    def is_charmm_HE1(atom, attached):
        ''' For alkene hydrogen: hydrogen of the type RHC=CR '''
        attached_sat = count_bonded(atom, lambda x: (x.elem == "C" and is_saturated(x)))
        attached_num_H = count_bonded(attached, lambda x: x.is_H)
        if not attached_sat and attached_num_H == 1: return True
        else: return False
    def is_charmm_HE2(atom, attached):
        ''' For alkene hydrogen: hydrogen of the type H2C=CR '''
        attached_sat = count_bonded(atom, lambda x: (x.elem == "C" and is_saturated(x)))
        attached_num_H = count_bonded(attached, lambda x: x.is_H)
        if not attached_sat and attached_num_H == 2: return True
        else: return False
    def is_charmm_HF1(atom, attached):
        ''' For Aliphatic H on fluorinated C: hydrogen attached to a carbon that is also attached to one fluorine'''
        if attached.elem == "C":
            attached_num_F = count_bonded(attached, lambda x: x.elem == "F")
            if attached_num_F == 1: return True
            else: return False
        else:return False
    def is_charmm_HF2(atom, attached):
        ''' For Aliphatic H on fluorinated C: hydrogen attached to a carbon that is also attached to two fluorine'''
        if attached.elem == "C":
            attached_num_F = count_bonded(attached, lambda x: x.elem == "F")
            if attached_num_F == 2: return True
            else: return False
        else:return False
    # ----- CARBON -----
    # ----- CARBON -----
    # ----- CARBON -----
    def is_charmm_CA(atom):
        ''' For aromatic carbon: carbon with 2 aromatic bonds '''
        num_hvy_bonds = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_aro_bonds = sum(1 for bond in atom.bonds if bond.order == Bond.AROMATIC)
        if num_aro_bonds == 2 and num_hvy_bonds >= 2: return True
        else: return False
    def is_charmm_CT(atom):
        ''' For aliphatic sp3 C without hydrogens: carbon with all single bonds none to hydrogen '''
        if is_saturated(atom) and len(atom.bonds) == 4:
            num_H = count_bonded(atom, lambda x: x.elem == "H")
            if num_H == 0: return True
            else: return False
        else: return False
    def is_charmm_CT1(atom):
        ''' For aliphatic sp3 C with 1 hydrogens: carbon with all single bonds one to hydrogen '''
        if is_saturated(atom) and len(atom.bonds) == 4:
            num_H = count_bonded(atom, lambda x: x.elem == "H")
            if num_H == 1: return True
            else: return False
        else: return False
    def is_charmm_CT2(atom):
        ''' For aliphatic sp3 C with 2 hydrogens: carbon with all single bonds two to hydrogen '''
        if is_saturated(atom) and len(atom.bonds) == 4:
            num_H = count_bonded(atom, lambda x: x.elem == "H")
            if num_H == 2: return True
            else: return False
        else: return False
    def is_charmm_CT3(atom):
        ''' For aliphatic sp3 C with 3 hydrogens: carbon with all single bonds three to hydrogen '''
        if is_saturated(atom) and len(atom.bonds) == 4:
            num_H = count_bonded(atom, lambda x: x.elem == "H")
            if num_H == 3: return True
            else: return False
        else: return False
    def is_charmm_CPT(atom):
        ''' For bridging carbons: carbon with 3 aromatic and 3 bonds to heavy atoms that bridges more than one aromatic ring'''
        num_hvy_bonds = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_aro_bonds = sum(1 for bond in atom.bonds if bond.order == Bond.AROMATIC)
        if num_aro_bonds == 3 and num_hvy_bonds == 3: return True
        else: return False
    def is_charmm_CS(atom):
        ''' For thiolate carbon: carbon attached to a "bare" sulfer ie. H3C*S(-) '''
        attached_S = [bond.a2 for bond in atom.bonds if bond.a2.elem == "S"]
        if len(attached_S) == 1:
            if len(attached_S[0].bonds) == 1: return True # does the sulfer only make one bond
            else: return False
        else: return False
    def is_charmm_CE1(atom):
        ''' For alkene carbon RHC*=CR: carbon double bonded to another carbond and single bonded to a hydorgen and something else '''
        attached_double = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.DOUBLE]
        num_H = count_bonded(atom, lambda x: x.is_H)
        if not is_saturated(atom) and len(attached_double) == 1 and num_H <= 1: return True
        else: return False
    def is_charmm_CE2(atom):
        ''' For alkene carbon H2C*=CR: carbon double bonded to another carbond and single bonded to two hydorgens '''
        attached_double = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.DOUBLE]
        num_H = count_bonded(atom, lambda x: x.is_H)
        if not is_saturated(atom) and len(attached_double) == 1 and num_H == 2: return True
    def is_charmm_CN(atom):
        ''' For carbon in cyano (nitrial) group: carbon triple bonded to a nitrogen '''
        num_trip_N = sum(1 for bond in atom.bonds if bond.order == Bond.TRIPLE and bond.a2.elem == "N")
        if num_trip_N == 1: return True
        else: return False
    def is_charmm_CPH1(atom):
        ''' For imidazole carbon: carbon in imidazol bonded to a N and a C ie. his CG and CD2 carbons '''
        if is_aromatic(atom):
            attached_N = [bond.a2 for bond in atom.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "N"]
            attached_C = [bond.a2 for bond in atom.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
            if len(attached_N) == 1 and len(attached_C) == 1:
                attached_C_N = [bond.a2 for bond in attached_C[0].bonds if bond.a2.elem == "N"]
                if len(attached_C_N) == 1: return True
                else: return False
            else: return False
        else: return False
    def is_charmm_CPH2(atom):
        ''' For imidazole carbon: carbon in imidiazole bonded to two N ie. his CE1 carbon '''
        if is_aromatic(atom):
            attached_N = [bond.a2 for bond in atom.bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N"]
            if len(attached_N) == 2: return True
            else: return False
        else: return False
    def is_charmm_CY(atom):
        ''' For indol carbon: carbon with aro bond to a bridge carbon (atom type CPT) '''
        attached_C = [bond.a2 for bond in atom.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
        for ac in attached_C:
            num_hvy_bonds = sum(1 for bond in ac.bonds if not bond.a2.elem == "H")
            num_aro_bonds = sum(1 for bond in ac.bonds if bond.order == Bond.AROMATIC)
            if num_aro_bonds == 3 and num_hvy_bonds == 3 and atom.ring_size == 5: return True
        return False
    def is_charmm_CC(atom):
        ''' For carboyl carbon: carbon double bonded to O and single bonded to an unprotonated O.
        Catches esters or carboxylic acids but not aldehyde, ketone'''
        attached_double = [bond.a2 for bond in atom.bonds if bond.a2.elem == "O" and bond.order == Bond.DOUBLE]
        attached_single = [bond.a2 for bond in atom.bonds if (bond.a2.elem == "O" or bond.a2.elem == "N") and bond.order == Bond.SINGLE ]
        if len(attached_double) == 1 and len(attached_single) == 1:
            if attached_single[0].elem == "O" and len(attached_single[0].bonds) == 1: return True
            elif attached_single[0].elem == "N" and len([True for bond in attached_single[0].bonds if bond.a2.elem == "H"]) == 2: return True
            else: return False
        else:return False
    def is_charmm_CD(atom):
        ''' For carboyl carbon: carbon double bonded to O and single bonded to a protonated O.
        Catches esters or carboxylic acids but not aldehyde, ketone, amides '''
        attached_double = [bond.a2 for bond in atom.bonds if bond.a2.elem == "O" and bond.order == Bond.DOUBLE]
        attached_single = [bond.a2 for bond in atom.bonds if bond.a2.elem == "O" and bond.order == Bond.SINGLE ]
        if len(attached_double) == 1 and len(attached_single) == 1:
            if len(attached_single[0].bonds) > 1: return True # O is bonded to something in addition to the carbon
            else: return False
        else:return False
    def is_charmm_C(atom):
        ''' For carboyl carbon: carbon double bonded to N or O and single bonded to two other atoms that are not O.
        Catches aldehyde, ketone, amides but not esters or carboxylic acids which are a special atom type '''
        attached_double = [bond.a2 for bond in atom.bonds if (bond.a2.elem == "O" or bond.a2.elem == "N") and bond.order == Bond.DOUBLE]
        attached_single = [bond.a2 for bond in atom.bonds if not bond.a2.elem == "O" and bond.order == Bond.SINGLE ]
        if len(attached_double) == 1 and len(attached_single) == 2: return True
        else:return False
    # ----- NITROGEN -----
    # ----- NITROGEN -----
    # ----- NITROGEN -----
    def is_charmm_N(atom):
        ''' For proline nitrogen: nitrogen with 3 heavy bonds '''
        num_hvy_bond = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_H = sum(1 for bond in atom.bonds if bond.a2.elem == "H")
        if num_hvy_bond == 3 and num_H == 0: return True
        else: return False
    def is_charmm_NR1(atom):
        ''' For neutral his protonated ring nitrogen: protonated aromatic nitrogen '''
        if is_aromatic(atom) and atom.ring_size == 5:
            attached_aro_C = [bond.a2 for bond in atom.bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
            attached_single = [bond.a2 for bond in atom.bonds if bond.order == Bond.SINGLE] # catches H and CH3
            if len(attached_aro_C) == 2 and len(attached_single) == 1:
                attached_aro_C_aro_N  = [bond.a2 for bond in attached_aro_C[0].bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N" and bond.a2.name != atom.name]
                attached_aro_C_aro_N += [bond.a2 for bond in attached_aro_C[1].bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N" and bond.a2.name != atom.name]
                if len(attached_aro_C_aro_N) == 1:
                    single = [bond.a2 for bond in attached_aro_C_aro_N[0].bonds if bond.order == Bond.SINGLE] # catches H and CH3
                    if len(single) == 0: return True
                    else: return False
                else: return False
        else:return False
    def is_charmm_NR2(atom):
        ''' For neutral his unprotonated ring nitrogen: unprotonated aromatic nitrogen '''
        if is_aromatic(atom) and atom.ring_size == 5:
            attached_aro_C = [bond.a2 for bond in atom.bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
            attached_single = [bond.a2 for bond in atom.bonds if bond.order == Bond.SINGLE] # catches H and CH3
            if len(attached_aro_C) == 2 and len(attached_single) == 0:
                attached_aro_C_aro_N  = [bond.a2 for bond in attached_aro_C[0].bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N" and bond.a2.name != atom.name]
                attached_aro_C_aro_N += [bond.a2 for bond in attached_aro_C[1].bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N" and bond.a2.name != atom.name]
                if len(attached_aro_C_aro_N) == 1:
                    single = [bond.a2 for bond in attached_aro_C_aro_N[0].bonds if bond.order == Bond.SINGLE] # catches H and CH3
                    if len(single) == 1: return True
                    else: return False
                else: return False
        else:return False
    def is_charmm_NR3(atom):
        ''' For charged his ring nitrogen: protonated nitrogen attached to a carbon attachec to another protonated nitrogen'''
        if is_aromatic(atom) and atom.ring_size == 5:
            attached_aro_C = [bond.a2 for bond in atom.bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
            attached_single = [bond.a2 for bond in atom.bonds if bond.order == Bond.SINGLE] # catches H and CH3
            if len(attached_aro_C) == 2 and len(attached_single) == 1:
                attached_aro_C_aro_N  = [bond.a2 for bond in attached_aro_C[0].bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N" and bond.a2.name != atom.name]
                attached_aro_C_aro_N += [bond.a2 for bond in attached_aro_C[1].bonds if bond.order == Bond.AROMATIC and bond.a2.elem == "N" and bond.a2.name != atom.name]
                if len(attached_aro_C_aro_N) == 1:
                    single = [bond.a2 for bond in attached_aro_C_aro_N[0].bonds if bond.order == Bond.SINGLE] # catches H and CH3
                    if len(single) == 1: return True
                    else: return False
                else: return False
        else:return False
    def is_charmm_NH1(atom):
        ''' For peptide nitrogen: nitrogen with 2 heavy bonds and 1 hydogen'''
        num_hvy_bond = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_H = sum(1 for bond in atom.bonds if bond.a2.elem == "H")
        if num_hvy_bond == 2 and num_H == 1: return True
        else: return False
    def is_charmm_NH2(atom):
        ''' For amide nitrogen: nitrogen with 1 heavy bond and 2 hydogen'''
        num_hvy_bond = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_H = sum(1 for bond in atom.bonds if bond.a2.elem == "H")
        if num_hvy_bond == 1 and num_H == 2: return True
        else: return False
    def is_charmm_NH3(atom):
        ''' For ammonium nitrogen: nitrogen with 1 heavy bond and 3 hydogen'''
        num_hvy_bond = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_H = sum(1 for bond in atom.bonds if bond.a2.elem == "H")
        if num_hvy_bond == 1 and num_H == 3: return True
        else: return False
    def is_charmm_NC2(atom):
        ''' For guanidinium nitrogen: nitrogen attached to carbon that is bonded to two other nitrogens'''
        attached_Cs = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C"]
        if len(attached_Cs) >= 1: #Some guanadinium N's are bound to multiple C's...
            for attached_C in attached_Cs:
                attached_C_N = [bond.a2 for bond in attached_C.bonds if bond.a2.elem == "N"]
                attached_C_a = [bond.a2 for bond in attached_C.bonds] #Gets rid of sp3 central carbons
                if len(attached_C_N) == 3 and len(attached_C_a) == 3: return True
            return False
        else: return False
    def is_charmm_NY(atom):
        ''' For indol nitrogen: nitrogen with aro bond to a bridge carbon (atom type CPT) '''
        attached_C = [bond.a2 for bond in atom.bonds if  bond.order == Bond.AROMATIC and bond.a2.elem == "C"]
        for ac in attached_C:
            num_hvy_bonds = sum(1 for bond in ac.bonds if not bond.a2.elem == "H")
            num_aro_bonds = sum(1 for bond in ac.bonds if bond.order == Bond.AROMATIC)
            if num_aro_bonds == 3 and num_hvy_bonds == 3 and atom.ring_size == 5: return True
        return False
    def is_charmm_NP(atom):
        ''' For nterm proline nitrogen: nitrogen with 2 heavy bond and 2 hydogen'''
        num_hvy_bond = sum(1 for bond in atom.bonds if not bond.a2.elem == "H")
        num_H = sum(1 for bond in atom.bonds if bond.a2.elem == "H")
        if num_hvy_bond == 2 and num_H == 2: return True
        else: return False
    def is_charmm_NC(atom):
        ''' For carbon in cyano (nitrial) group: carbon triple bonded to a nitrogen '''
        num_trip_C = sum(1 for bond in atom.bonds if bond.order == Bond.TRIPLE and bond.a2.elem == "C")
        if num_trip_C == 1: return True
        else: return False
    # ----- OXYGEN -----
    # ----- OXYGEN -----
    # ----- OXYGEN -----
    def is_charmm_O(atom):
        ''' For carbonyl oxygen: oxygen only double bonded to a carbon that is not bonded to another oxygen '''
        attached_double_C = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.DOUBLE]
        attached_single_C = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.SINGLE]
        if len(attached_double_C) == 1 and len(attached_single_C) == 0:
            attached_double_C_O = [bond.a2 for bond in attached_double_C[0].bonds if bond.a2.elem == "O" and bond.a2.name != atom.name]
            if len(attached_double_C_O) == 0: return True
            else: return False
        else: return False
    def is_charmm_OB(atom):
        ''' For carbonyl oxygen in acetic acid: oxygen only double bonded to a carbon that is single bonded to another oxygen that itself has a bond to something else'''
        attached_double_C = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.DOUBLE]
        attached_single_C = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.SINGLE]
        if len(attached_double_C) == 1 and len(attached_single_C) == 0:
            attached_double_C_O = [bond.a2 for bond in attached_double_C[0].bonds if bond.a2.elem == "O" and bond.a2.name != atom.name]
            if len(attached_double_C_O) == 1:
                attached_double_C_O_any = [bond.a2 for bond in attached_double_C_O[0].bonds if bond.a2.name != attached_double_C[0].name]
                if len(attached_double_C_O_any) > 0: return True
                else: return False
            else: return False
        else: return False
    def is_charmm_OC(atom):
        ''' For carboxylate oxygen: oxygen only (A)single/(B)double bonded to a carbon that is (A)double/(B)single bonded to another bare oxygen that is '''
        attached_double_C = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.DOUBLE]
        attached_single_C = [bond.a2 for bond in atom.bonds if bond.a2.elem == "C" and bond.order == Bond.SINGLE]
        # case A
        if len(attached_double_C) == 1 and len(attached_single_C) == 0:
            attached_double_C_O = [bond.a2 for bond in attached_double_C[0].bonds if bond.a2.elem == "O" and bond.a2.name != atom.name]
            if len(attached_double_C_O) == 1:
                attached_double_C_O_any = [bond.a2 for bond in attached_double_C_O[0].bonds if bond.a2.name != attached_double_C[0].name]
                if len(attached_double_C_O_any) == 0: return True
                else: return False
            else: return False
        # case B
        elif len(attached_double_C) == 0 and len(attached_single_C) == 1:
            attached_single_C_O = [bond.a2 for bond in attached_single_C[0].bonds if bond.a2.elem == "O" and bond.a2.name != atom.name]
            if len(attached_single_C_O) == 1:
                attached_single_C_O_any = [bond.a2 for bond in attached_single_C_O[0].bonds if bond.a2.name != attached_single_C[0].name]
                if len(attached_single_C_O_any) == 0: return True
                else: return False
            else: return False
        else: return False
    def is_charmm_OH1(atom):
        ''' For hydroxyl oxygen: oxygen bonded to a hydrogen and somthing else '''
        attached_H = [bond.a2 for bond in atom.bonds if bond.a2.elem == "H"]
        attached_other = [bond.a2 for bond in atom.bonds if not bond.a2.elem == "H"]
        if len(attached_H) == 1 and len(attached_other) >= 1: return True
        else: return False
    def is_charmm_OS(atom):
        ''' For ester oxygen: oxygen with 2 heavy atom single bonds, think this may work for ether too'''
        attached_H = [bond.a2 for bond in atom.bonds if bond.a2.elem == "H"]
        attached_other = [bond.a2 for bond in atom.bonds if not bond.a2.elem == "H"]
        if len(attached_H) == 0 and len(attached_other) >= 2: return True
        else: return False
    # ----- SULFER -----
    # ----- SULFER -----
    # ----- SULFER -----
    def is_charmm_S(atom):
        ''' For sulfur: sulfur with more than one bond '''
        num_bonds = len(atom.bonds)
        if num_bonds > 1: return True
        else: return False
    def is_charmm_SS(atom):
        ''' For thiolate sulfur: bare sulfur connected to a carbon '''
        num_bonds = len(atom.bonds)
        if num_bonds == 1: return True
        else: return False
    # main loop
    for a in atoms:
        if a.is_virtual: a.mm_type = "VIRT"
        elif a.elem == "H":
            assert(len(a.bonds) == 1) # hydrogen should only be attached to one atom
            at = a.bonds[0].a2
            #HC comes first because of delocalized/aromatic NH2
            if   is_charmm_HC(a, at):  a.mm_type = "HC" 
            elif is_charmm_H(a, at):   a.mm_type = "H"
            elif is_charmm_HS(a, at):  a.mm_type = "HS"
            elif is_charmm_HB(a, at, peptoid):  a.mm_type = "HB"
            elif is_charmm_HA(a, at):  a.mm_type = "HA"
            elif is_charmm_HR1(a, at): a.mm_type = "HR1"
            elif is_charmm_HR2(a, at): a.mm_type = "HR2"
            elif is_charmm_HR3(a, at): a.mm_type = "HR3"
            elif is_charmm_HP(a, at):  a.mm_type = "HP"
            elif is_charmm_HE1(a, at): a.mm_type = "HE1"
            elif is_charmm_HE2(a, at): a.mm_type = "HE2"
            elif is_charmm_HF1(a, at): a.mm_type = "HF1"
            elif is_charmm_HF2(a, at): a.mm_type = "HF2"
            else: a.mm_type = "X"
        elif a.elem == "C" :
            if   is_charmm_CT(a):   a.mm_type = "CT1"
            elif is_charmm_CT1(a):  a.mm_type = "CT1"
            elif is_charmm_CT2(a):  a.mm_type = "CT2"
            elif is_charmm_CT3(a):  a.mm_type = "CT3"
            elif is_charmm_CPH1(a): a.mm_type = "CPH1"
            elif is_charmm_CPH2(a): a.mm_type = "CPH2"
            elif is_charmm_CPT(a):  a.mm_type = "CPT"
            elif is_charmm_CY(a):   a.mm_type = "CY"
            elif is_charmm_CA(a):   a.mm_type = "CA"
            elif is_charmm_CC(a):   a.mm_type = "CC"
            elif is_charmm_CD(a):   a.mm_type = "CD"
            elif is_charmm_CS(a):   a.mm_type = "CS"
            elif is_charmm_CE1(a):  a.mm_type = "CE1"
            elif is_charmm_CE2(a):  a.mm_type = "CE2"
            elif is_charmm_CN(a):   a.mm_type = "CN"
            elif is_charmm_C(a):    a.mm_type = "C"
            else: a.mm_type = "X"
        elif a.elem == "N" :
            if   is_charmm_NY(a):  a.mm_type = "NY"
            elif is_charmm_NR1(a): a.mm_type = "NR1"
            elif is_charmm_NR2(a): a.mm_type = "NR2"
            elif is_charmm_NR3(a): a.mm_type = "NR3"
            elif is_charmm_N(a):   a.mm_type = "N"
            elif is_charmm_NC2(a): a.mm_type = "NC2"
            elif is_charmm_NP(a):  a.mm_type = "NP"
            elif is_charmm_NC(a):  a.mm_type = "NC"
            #These are so generically defined they'll override more specific N types
            elif is_charmm_NH1(a): a.mm_type = "NH1"
            elif is_charmm_NH2(a): a.mm_type = "NH2"
            elif is_charmm_NH3(a): a.mm_type = "NH3"
            else: a.mm_type = "X"
        elif a.elem == "O" :
            if   is_charmm_O(a):   a.mm_type = "O"
            elif is_charmm_OB(a):  a.mm_type = "OB"
            elif is_charmm_OC(a):  a.mm_type = "OC"
            elif is_charmm_OH1(a): a.mm_type = "OH1"
            elif is_charmm_OS(a):  a.mm_type = "OS"
            else: a.mm_type = "X"
        elif a.elem == "S" :
            if   is_charmm_SS(a):  a.mm_type = "SS"
            elif is_charmm_S(a):   a.mm_type = "S"
            else: a.mm_type = "X"
        elif a.elem == "P" : a.mm_type = "P"
        elif a.elem == "F" : a.mm_type = "F1"
        elif a.elem == "CL": a.mm_type = "CL"
        elif a.elem == "BR": a.mm_type = "BR"
        elif a.elem == "I" : a.mm_type = "I"
        elif a.elem == "NA": a.mm_type = "NA"
        elif a.elem == "K" : a.mm_type = "K"
        elif a.elem == "MG": a.mm_type = "MG"
        elif a.elem == "FE": a.mm_type = "FE"
        elif a.elem == "CA": a.mm_type = "CA"
        elif a.elem == "ZN": a.mm_type = "ZN"
        else: a.mm_type = " X  " # this at least seems to be a legal MM atom type

def read_parital_charge_input(input_fname):
    f = open(input_fname, "r")
    partial_charges = []
    for x in f:
        partial_charges.append(x.split())
    return partial_charges

def assign_partial_charges_from_values(molfile, partial_charges, net_charge=0):
    atoms = molfile.atoms
    # if the input partial charge file is provided, use the information in that file
    assert len(atoms) == len(partial_charges), "Error: the numbers of partial charge values and atom are not the same"
    current_net_charge = 0.0
    for charge in partial_charges:
        current_net_charge += float(charge[2])
        #print( "charge is %s" % current_net_charge )
    assert abs(net_charge - current_net_charge) < 1e-4, "Error: sum of partial charges %1.3f is not equal to the net charge %1.3f " % (current_net_charge, net_charge)
    for i in range(len(atoms)):
        assert atoms[i].elem == partial_charges[i][1].upper(),"Error: the elements %s in the sdf file and %s the partial charge file doesn't match" % (atoms[i].elem, partial_charges[i][1].upper())
        atoms[i].partial_charge = float(partial_charges[i][2])

def assign_partial_charges(atoms, partial_charges, net_charge=0.0):
    '''Assigns Rosetta standard partial charges, then
    corrects them so they sum to the desired net charge.
    Correction is distributed equally among all atoms.

    If non-zero partial charges were already assigned, no change is made.
    '''
    # If the partial charges file is not provided
    if partial_charges == None:
        
        std_charges = { # from Rosetta++ aaproperties_pack.cc
            "CNH2" : 0.550,
            "COO " : 0.620,
            "CH0 " : 0.075,
            "CH1 " : -0.090,
            "CH2 " : -0.180,
            "CH3 " : -0.270,
            "aroC" : -0.115,
            "Ntrp" : -0.610,
            "Nhis" : -0.530,
            "NH2O" : -0.470,
            "Nlys" : -0.620,
            "Narg" : -0.750,
            "Npro" : -0.370,
            "OH  " : -0.660,
            "Oaro" : -0.660, # copied from OH
            "ONH2" : -0.550,
            "OOC " : -0.760,
            "S   " : -0.090,
            "SH1 " : -0.230,
            "Nbb " : -0.470,
            "CAbb" : 0.070,
            "CObb" : 0.510,
            "OCbb" : -0.510,
            "Phos" : 1.500,
            "Hpol" : 0.430,
            "Hapo" : 0.095,
            "Haro" : 0.115,
            "HNbb" : 0.310,
            "H2O " : 0.000,
            "HS  " : 0.160,
            "F   " : -0.250,
            "Cl  " : -0.130,
            "Br  " : -0.100,
            "I   " : -0.090,
            "Zn2p" : 2.000,
            "Fe2p" : 2.000,
            "Fe3p" : 3.000,
            "Mg2p" : 2.000,
            "Ca2p" : 2.000,
            "Na1p" : 1.000,
            "K1p " : 1.000,
            "Bsp2" : 0.020,
            "VIRT" : 0.000,
        }

        for a in atoms:
            a.partial_charge = std_charges[ a.ros_type ]
    null_charge = [a for a in atoms if a.partial_charge is None]
    if 0 < len(null_charge): 
        if len(null_charge) < len(atoms):
            raise ValueError("Only some partial charges were assigned -- must be all or none.")
        else:
            return
    # We only want to operate on non-virtual atoms now:
    atoms = [a for a in atoms if not (a.poly_ignore or a.is_virtual or a.poly_lower or a.poly_upper)]
    curr_net_charge = sum(a.partial_charge for a in atoms)
    # check if the current sum of all partial charges is the same with the net charge
    charge_correction = (net_charge - curr_net_charge) / len(atoms)
    if abs(charge_correction) > 1e-4:
        print("Total naive charge %.3f, desired charge %.3f, offsetting all atoms by %.3f" %
              (curr_net_charge, net_charge, charge_correction)
          )
        curr_net_charge = 0.0
        for a in atoms:
            a.partial_charge += charge_correction
            curr_net_charge += a.partial_charge
        assert abs(net_charge - curr_net_charge) < 1e-4, "charge correction failed"

def compare_molfiles(m1, m2):
    '''
    Check if two molfiles represent the same molecule.
    Check if the molecule size and the atom order match
    '''
    # TODO: raise exception instead of assert as we can skip the molfile if things go wrong
    if len(m1.atoms) != len(m2.atoms):
        print("Two molecules have different number of atoms %d %d" % (len(m1.atoms), len(m2.atoms)))
        return False
    if len(m1.bonds) != len(m2.bonds):
        print("Two molecules have different number of bonds")
        return False
    for inx, atom in enumerate(m1.atoms):
        if atom.elem != m2.atoms[inx].elem:
            print("Two molecules have unmatched atoms %d %s %s" % (inx+1, atom.elem, m2.atoms[inx].elem))
            return False
    return True

def copy_atom_and_bond_info(original, copy):
    '''
    Copies atom and bond information between from original to copy molefiles
    '''
    assert compare_molfiles(original, copy), \
    "Error: 2 molecules are different, expect the same molecule"
    #add_fields_to_atoms(copy.atoms)
    for inx, atom in enumerate(original.atoms):
        catom = copy.atoms[inx]
        catom.orig_name = atom.orig_name # for kinemage output
        catom.pdb_name = atom.pdb_name       # PDB style atom name
        catom.ros_type = atom.ros_type         # Rosetta atom type
        catom.mm_type =  atom.mm_type        # Molec. mechan. atom type
        catom.is_virtual = atom.is_virtual   # for atom typing and charge assignment
        catom.rigid_id = atom.rigid_id          # non-zero id for regions with no rotatable bonds; may span fragments
        catom.fragment_id = atom.fragment_id       # non-zero id for fragments after bond breaking
        #catom.conn_bonds = []       # list of cross-fragment bonds to this atom
        catom.is_root = atom.is_root       # for atom tree
        #catom.parent =  atom.parent        # for atom tree
        #catom.children = []         # for atom tree
        #catom.stub1 = None          # for internal coords, derived from atom tree
        #catom.stub2 = None          # for internal coords, derived from atom tree
        #catom.stub3 = None          # for internal coords, derived from atom tree
        #catom.input_stub1 = None    # for internal coords, derived from atom tree
        #catom.input_stub2 = None    # for internal coords, derived from atom tree
        #catom.input_stub3 = None    # for internal coords, derived from atom tree
        #catom.d = 0.0               # distance to input_stub1
        #catom.theta = 0.0           # 180 - angle with input_stub2 (degrees)
        #catom.phi = 0.0             # dihedral from input_stub3 (degrees)
        catom.poly_upper = atom.poly_upper    # is upper connect atom for polymer residue type
        catom.poly_lower = atom.poly_lower    # is lower connect atom for polymer residue type
        catom.poly_n_bb = atom.poly_n_bb     # is backbone nitrogen for polymer residue type
        catom.poly_ca_bb = atom.poly_ca_bb    # is backbone alpha carbon for polymer residue type
        catom.poly_c_bb = atom.poly_c_bb     # is backbone carbonyl carbon for polymer residue type
        catom.poly_o_bb = atom.poly_o_bb     # is backbone carbonyl oxygen for polymer residue type
        catom.poly_backbone = atom.poly_backbone # convience boolean
        catom.poly_ignore = atom.poly_ignore   # convience boolean
        catom.pdb_prefix_num = atom.pdb_prefix_num
        catom.pdb_elem = catom.elem  # blah
        catom.pdb_greek_dist = atom.pdb_greek_dist  # bond distance?
        catom.pdb_postfix_num = atom.pdb_postfix_num
        # copy bond info over
        #add_fields_to_bonds(copy.bonds)
        for inx, bond in enumerate(original.bonds):
            cbond = copy.bonds[inx]
            cbond.can_rotate = bond.can_rotate     # true for single bonds not in rings
            cbond.is_proton_chi = bond.is_proton_chi  # true for bonds that rotate OH hydrogens, etc
            cbond.connection_id = bond.connection_id      # non-zero id if bond crosses fragments
            # Remember we have to update mirror too!
            cbond.mirror.can_rotate      = bond.can_rotate
            cbond.mirror.is_proton_chi   = bond.is_proton_chi
            cbond.mirror.connection_id   = bond.connection_id
            cbond.poly_ignore = bond.poly_ignore  # convience boolean
        # after atoms and bonds info are copied over,
        # we can now reorder the atoms of both original and copy
