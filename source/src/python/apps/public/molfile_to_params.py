#!/usr/bin/env python
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
'''
Functions and executable for taking a ligand from an MDL Molfile
and splitting it into one or more .params files for Minirosetta.
See main() for usage or run with --help.

Author: Ian W. Davis
'''
import os, sys, copy, random #{{{
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

from optparse import OptionParser, IndentedHelpFormatter

try: set
except: from sets import Set as set

# Magic spell to make sure Rosetta python libs are on the PYTHONPATH:
sys.path.append( os.path.dirname(os.path.dirname( os.path.abspath(sys.path[0]) )) )

from rosetta_py.io.mdl_molfile import *
from rosetta_py.utility.rankorder import argmin, order
from rosetta_py.utility import r3

# Features from Python 2.5 that we want to use:
try: any
except:
    def any(itr):
        for el in itr:
            if el: return True
        return False

try: all
except:
    def all(itr):
        for el in itr:
            if not el: return False
        return True

# Better handle multiple paragraph descriptions.
class PreformattedDescFormatter (IndentedHelpFormatter):
    def format_description(self, description):
        return description.strip() + "\n" # Remove leading/trailing whitespace

# Quick and dirty memoize to speed up neighbor atom calculations
def memoize(f):
    cache = {}
    def decorated(*args):
        if args not in cache:
            cache[args] = f(*args)
        return cache[args]
    return decorated
#}}}

def mark_fragments(molfile): #{{{
    # Part 0:  initialize a few specific atom and bond properties
    for atom in molfile.atoms:
        atom.is_nbr = False     # did user input request this as nbr_atom?
    for bond in molfile.bonds:
        bond.break_me = False   # did user input request to break this bond?
        bond.mirror.break_me = bond.break_me
    # Part 1:  find bonds explicitly marked to break
    for line in molfile.footer:
        if not line.startswith("M SPLT"): continue
        fields = line.split()
        atom1 = molfile.atoms[int(fields[2]) - 1]
        atom2 = molfile.atoms[int(fields[3]) - 1]
        for b in molfile.bonds:
            if((b.a1 == atom1 and b.a2 == atom2)
            or (b.a1 == atom2 and b.a2 == atom1)):
                b.break_me = True
                b.mirror.break_me = b.break_me
                break
        else:
            raise ValueError("Cannot find bond to split between %s and %s" % (atom1.name, atom2.name))
    # Part 2:  break all bonds around a "fragment" atom and create virtual atoms
    ctr_lines = [ line for line in molfile.footer if line.startswith("M FRAG") ]
    if len(ctr_lines) == 0:
        return
    elif len(ctr_lines) > 1:
        raise ValueError("Cannot fragment ligand around more than one atom")
    ctr_atom = molfile.atoms[int(ctr_lines[0].split()[2]) - 1]
    if ctr_atom.is_H:
        raise ValueError("Cannot fragment ligand around a hydrogen atom")
    elif len(ctr_atom.heavy_bonds) < 2:
        raise ValueError("Cannot fragment ligand around an atom with less than 2 heavy-atom neighbors")
    ctr_atom.is_root = True # not sure this is necessary
    ctr_atom.is_nbr = True  # the virtualized center atom must be the nbr_atom for ALL fragments
    def add_virtual(a1, a2, order):
        v2 = a2.copy()
        v2.name = "V" # for "virtual"
        #v2.sybyl_type = "Du" # for "Dummy atom"
        # Keeps dummies from being on top of each other and/or in straight lines,
        # which would lead to NaNs in the atom tree inside Mini.
        r3.add( v2, r3.Triple(random.uniform(-0.01, 0.01), random.uniform(-0.01, 0.01), random.uniform(-0.01, 0.01)), v2 )
        molfile.atoms.append(v2)
        b = Bond(a1, v2, order) # automatically updates bonds[] and heavy_bonds[] for the atoms
        b.break_me = False
        b.mirror.break_me = b.break_me
        molfile.bonds.append(b)
        return v2
    # ctr_atom becomes an island, bonded to virtual copies of its original neighbors
    # each neighbor gets a virtual copy of ctr_atom and all ctr_atom's other neighbors
    ctr_bonds = list(ctr_atom.heavy_bonds) # iterate over a copy, b/c it will be modified in the loop
    for bond in ctr_bonds:
        bond.break_me = True
        bond.mirror.break_me = bond.break_me
        assert bond.a1 is ctr_atom
        add_virtual(bond.a1, bond.a2, bond.order)
        vctr = add_virtual(bond.a2, bond.a1, bond.order)
        vctr.is_root = True
        vctr.is_nbr = True
        for b in ctr_bonds:
            if b is bond: continue
            add_virtual(vctr, b.a2, b.order)
#}}}
def add_fields_to_atoms(atoms): #{{{
    '''Adds a bunch of member variable placeholders that we use.'''
    for atom in atoms:
        atom.orig_name = atom.name # for kinemage output
        atom.ros_type = ""      # Rosetta atom type
        atom.mm_type = ""       # Molec. mechan. atom type
        atom.cen_type = ""      # Rosetta centroid atom type
        atom.is_virtual = False # for atom typing and charge assignment
        atom.rigid_id = 0       # non-zero id for regions with no rotatable bonds; may span fragments
        atom.fragment_id = 0    # non-zero id for fragments after bond breaking
        atom.conn_bonds = []    # list of cross-fragment bonds to this atom
        atom.is_root = False    # for atom tree
        atom.parent = None      # for atom tree
        atom.children = []      # for atom tree
        atom.stub1 = None       # for internal coords, derived from atom tree
        atom.stub2 = None       # for internal coords, derived from atom tree
        atom.stub3 = None       # for internal coords, derived from atom tree
        atom.input_stub1 = None # for internal coords, derived from atom tree
        atom.input_stub2 = None # for internal coords, derived from atom tree
        atom.input_stub3 = None # for internal coords, derived from atom tree
        atom.d = 0.0            # distance to input_stub1
        atom.theta = 0.0        # 180 - angle with input_stub2 (degrees)
        atom.phi = 0.0          # dihedral from input_stub3 (degrees)
#}}}
def add_fields_to_bonds(bonds): #{{{
    '''Adds a bunch of member variable placeholders that we use.'''
    for bond in bonds:
        bond.can_rotate = False     # true for single bonds not in rings
        bond.is_proton_chi = False  # true for bonds that rotate OH hydrogens, etc
        bond.connection_id = 0      # non-zero id if bond crosses fragments
        # Remember we have to update mirror too!
        bond.mirror.can_rotate      = bond.can_rotate
        bond.mirror.is_proton_chi   = bond.is_proton_chi
        bond.mirror.connection_id   = bond.connection_id
#}}}
def find_virtual_atoms(atoms): #{{{
    '''Atoms whose names start with "V" are virtual, used in enzyme design, etc.'''
    for atom in atoms:
        if atom.name.startswith("V") or atom.name.startswith("X"):
            atom.is_virtual = True
            atom.elem = "X" # "V" is a real element (Vanadium)
#}}}
def check_bond_count(atoms): #{{{
    '''Safety check for structures with stupid numbers of bonds to atoms, at Florian's request.'''
    def valence(bond):
        if bond.order == Bond.AROMATIC: return 1.5
        else: return float(bond.order)
    for atom in atoms:
        if atom.is_H and len(atom.bonds) > 1:
            raise ValueError("Atom %s is a hydrogen with >1 bonds" % atom.name)
        # Valence 4.5 for e.g. carbon comes up in peptide bonds and at the joins in multi-ring systems.
        if sum([valence(bond) for bond in atom.bonds]) > 4.5:
            print "WARNING:  atom %s has valence > 4" % atom.name
#}}}
def check_aromaticity(bonds): #{{{
    '''Safety check for Kekule structures (alternating single/double bonds)
    rather than bonds described explicitly as aromatic.'''
    aro_bonds = len([b for b in bonds if b.order == Bond.AROMATIC])
    dbl_bonds = len([b for b in bonds if b.order == Bond.DOUBLE])
    if aro_bonds == 0 and dbl_bonds > 0:
        print "WARNING:  structure contains double bonds but no aromatic bonds"
        print "  Aromatic bonds must be identified explicitly --"
        print "  alternating single/double bonds (Kekule structure) won't cut it."
        print "  This warning does not apply to you if your molecule really isn't aromatic."
#}}}
def check_hydrogens(atoms): #{{{
    '''Safety check for from PDB structures.
    If you convert an X-ray structure of a ligand to a mol file, chances are it won't have hydrogens added.
    This will probably mess up the parameters.'''
    for a in atoms:
        if "H" == a.elem:
            return   # At least one hydrogen - we're probably fine.
    print "WARNING:  structure does not contain any hydrogens"
    print "  Hydrogens aren't automatically added. --"
    print "  Check your PDB -> mol conversion program for hydrogen-addition options."
    print "  This warning does not apply to you if your molecule shouldn't contain any hydrogens."
#}}}
def assign_rosetta_types(atoms): #{{{
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
        elif "H"  == a.elem:
            num_aro_C = count_bonded(a, lambda x: (x.elem == "C" and is_aromatic(x)) or x.ros_type == "aroC")
            num_NOS = count_bonded(a, lambda x: x.elem == "N" or x.elem == "O" or x.elem == "S")
            if num_NOS >= 1:    a.ros_type = "Hpol"
            elif num_aro_C >=1: a.ros_type = "Haro"
            else:               a.ros_type = "Hapo"
        elif "C"  == a.elem:
            if is_saturated(a):
                num_H = count_bonded(a, lambda x: x.is_H)
                if num_H >= 3:      a.ros_type = "CH3 "
                elif num_H == 2:    a.ros_type = "CH2 "
                else:               a.ros_type = "CH1 "
            else:
                num_dbl_nonO = 0; num_aro_nonO = 0; num_aro_N = 0;
                a_bonds = [b for b in a.bonds if not b.a2.is_virtual]
                for bond in a_bonds:
                    if bond.order == Bond.DOUBLE:
                        if bond.a2.elem != "O": num_dbl_nonO += 1
                    elif bond.order == Bond.AROMATIC:
                        if bond.a2.elem != "O": num_aro_nonO += 1
                        if bond.a2.elem == "N": num_aro_N += 1 # really if, not elif
                #print i+1, a.name, num_aro_nonO, num_dbl_nonO, num_aro_N
                if num_aro_nonO >= 2:   a.ros_type = "aroC"
                elif num_dbl_nonO >= 1: a.ros_type = "aroC"
                elif num_aro_N >= 1:    a.ros_type = "CNH2"
                else:                   a.ros_type = "COO "
        elif "N"  == a.elem:
            num_H = count_bonded(a, lambda x: x.is_H)
            heavy_nbrs = count_bonded(a, lambda x: not x.is_H)
            #assert num_H + heavy_nbrs == len(a.bonds) -- this fails if there are virtual atoms involved
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
        elif "S"  == a.elem: a.ros_type = "S   "
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
        else: raise ValueError("Unknown element '%s'" % a.elem)
#}}}
def assign_mm_types(atoms,mm_as_virt): #{{{
    '''For now, just fills in dummy values.'''
    for a in atoms:
        if a.is_virtual or mm_as_virt:
            a.mm_type = "VIRT"
        else:
            a.mm_type = " X  " # this at least seems to be a legal MM atom type
#}}}
def assign_centroid_types(atoms): #{{{
    '''Uses Nbb for donors, OCbb for acceptors, and CAbb for everything else (but H)'''
    Hbb = ["Hpol", "HNbb"] # polar H must be present for H-bonding?
    Nbb = ["Ntrp", "NH2O", "Nlys", "Narg", "Npro", "Nbb "]
    Obb = ["Nhis", "OH  ", "ONH2", "OOC ", "Oaro", "OCbb"]
    for a in atoms:
        if a.ros_type in Hbb: a.cen_type = "HNbb"
        elif a.ros_type in Nbb: a.cen_type = "Nbb "
        elif a.ros_type in Obb: a.cen_type = "OCbb"
        elif a.is_H: a.cen_type = None # nonpolar hydrogens get ignored
        else: a.cen_type = "CAbb"
#}}}
def setup_amino_acid(atoms, molfiles): #{{{
    '''Tries to reconfigure the input as a (modified) amino acid type,
    following the implicit Rosetta rules for these residues.
    '''
    # Sort atoms so some are first, some are last, and rest are in the middle (in input order).
    # Fortunately Python's sort() is stable.
    atom_prio = { " N  ":-6, " CA ":-5, " C  ":-4, " O  ":-2, " CB ":-1, #, "UPPER":-3
        # ... default atoms have priority 0 ...
        " HA ":1, " H  ":3 } #, "LOWER":2
    # Have to sort all molfiles the same way to maintain correspondances
    #atoms.sort( key=lambda a: atom_prio.get(a.name, 0) )
    atom_order = order([atom_prio.get(a.name, 0) for a in atoms])
    for molfile in molfiles:
        molfile.atoms = [ molfile.atoms[i] for i in atom_order ]
    # Within each atom, sort bonds in same order
    for a in atoms:
        a.bonds.sort( key=lambda b: atom_prio.get(b.a2.name, 0) )
    # These two things should (?) be enough to make the .params entries come out in the right order
    # Now, check to make sure all the required atoms are present
    atom_names = dict( (a.name, a) for a in atoms )
    if not all( name in atom_names for name in atom_prio.iterkeys() ):
        missing = set(atom_prio.iterkeys()) - set(atom_names.iterkeys())
        raise ValueError("Input is missing essential atom(s) for amino acids: %s" % ", ".join(missing))
    # Make N the root
    atom_names[" N  "].is_root = True
    # Set atom types (Rosetta and MM)
    atom_types = {
        " N  ":  ("N", "Nbb ", "NH1 "),
        " CA ": ("C", "CAbb", "CT1 "),
        " C  ":  ("C", "CObb", "C   "),
        " O  ":  ("O", "OCbb", "O   "),
        " H  ":  ("H", "HNbb", "H   "),
        " HA ": ("H", "Hapo", "HB  "),
    }
    for name, atype in atom_types.iteritems():
        atom_names[name].elem       = atype[0]
        atom_names[name].ros_type   = atype[1]
        atom_names[name].mm_type    = atype[2]
#}}}
def assign_partial_charges(atoms, net_charge=0.0, recharge=False): #{{{
    '''Assigns Rosetta standard partial charges, then
    corrects them so they sum to the desired net charge.
    Correction is distributed equally among all atoms.

    If non-zero partial charges were already assigned, no change is made.
    '''
    if not recharge:
        null_charge = [a for a in atoms if a.partial_charge is None]
        abs_charge = sum(abs(a.partial_charge) for a in atoms if a.partial_charge is not None)
        if len(null_charge) == 0 and abs_charge > 0:
            net_charge = sum(a.partial_charge for a in atoms if a.partial_charge is not None)
            print "Partial charges already fully assigned, no changes made; net charge %.3f" % net_charge
            return
        elif 0 < len(null_charge) and len(null_charge) < len(atoms):
            raise ValueError("Only some partial charges were assigned -- must be all or none.")
    else:
        print "Ignoring partial charges in file; net charge set to %.3f" % net_charge
    std_charges = { # from Rosetta++ aaproperties_pack.cc
        "CNH2" : 0.550,
        "COO " : 0.620,
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
        "S   " : -0.160,
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
        "VIRT" : 0.000,
    }
    curr_net_charge = 0.0
    for a in atoms:
        a.partial_charge = std_charges[ a.ros_type ]
        curr_net_charge += a.partial_charge
    # We only want to operate on non-virtual atoms now:
    atoms = [a for a in atoms if not a.is_virtual]
    charge_correction = (net_charge - curr_net_charge) / len(atoms)
    print "Total naive charge %.3f, desired charge %.3f, offsetting all atoms by %.3f" % (curr_net_charge, net_charge, charge_correction)
    curr_net_charge = 0.0
    for a in atoms:
        a.partial_charge += charge_correction
        curr_net_charge += a.partial_charge
    assert( abs(net_charge - curr_net_charge) < 1e-3 )
#}}}
def assign_rotatable_bonds(bonds): #{{{
    '''Rotatable bonds are single bonds outside of rings
    with heavy atoms attached to both ends (i.e. not methyls) or non-C with one H.'''
    def is_ok(a):
        bonds_to_heavy = len(a.heavy_bonds)
        bonds_to_hydro = len([b for b in a.bonds if b.a2.is_H])
        return bonds_to_heavy >= 2 or (
            bonds_to_heavy == 1 and a.heavy_bonds[0].order == Bond.SINGLE
            and bonds_to_hydro == 1 and a.elem != "C"
        )
    for b in bonds:
        b.can_rotate = (
            b.order == Bond.SINGLE
            and not b.is_ring
            and is_ok(b.a1)
            and is_ok(b.a2)
        )
        b.is_proton_chi = (
            b.can_rotate
            and (len(b.a1.heavy_bonds) == 1 or len(b.a2.heavy_bonds) == 1)
        )
        b.mirror.can_rotate     = b.can_rotate
        b.mirror.is_proton_chi  = b.is_proton_chi
#}}}
def assign_rigid_ids(atoms): #{{{
    '''Groups atoms that are connected in rigid units, i.e. no rotatable bonds.'''
    # Iterate through atoms, assigning them to rigids via depth first search
    def assign_to_rigid(atom, rig_id):
        if atom.rigid_id == rig_id: return # already visited this one
        assert(atom.rigid_id == 0)
        atom.rigid_id = rig_id
        # Recurse over all non-rotatable bonds and assign those atoms to same rigid
        for bond in atom.bonds:
            if not bond.can_rotate:
                assign_to_rigid(bond.a2, rig_id)
    num_rig_id = 0
    for atom in atoms:
        # We'll skip all but the first atom in each rigid as all atoms
        # in the rigid will be assigned by assign_to_rigid()
        if atom.rigid_id == 0:
            num_rig_id += 1
            assign_to_rigid(atom, num_rig_id)
#}}}
def fragment_ligand(molfile): #{{{
    '''Sets Atom.fragment_id, Atom.conn_bonds, and Bond.connection_id.
    Returns number of fragments created, i.e. the largest valid fragment id.'''
    remaining_bonds = set(molfile.bonds) # copy
    # Delete all split bonds from remaining_bonds
    # and number the connections they leave behind
    num_conn_id = 0
    for b in molfile.bonds:
        if not b.break_me: continue
        if b.can_rotate:
            print "WARNING: spliting ROTATABLE bond between %s and %s" % (b.a1.name, b.a2.name)
        print "Split bond between %s and %s" % (b.a1.name, b.a2.name)
        num_conn_id += 1
        b.connection_id = num_conn_id
        b.mirror.connection_id = b.connection_id
        b.a1.conn_bonds.append(b)
        b.a2.conn_bonds.append(b.mirror)
        remaining_bonds.remove(b)
    # Even though no single fragment can have more than 9 connections (CONN1 - CONN9),
    # the set of all fragments together may have more than 9 connections.
    #
    # Iterate through atoms, assigning them to fragments via depth first search
    def assign_to_fragment(atom, frag_id):
        if atom.fragment_id == frag_id: return # already visited this one
        assert(atom.fragment_id == 0)
        atom.fragment_id = frag_id
        # Recurse over all unbroken bonds and assign those atoms to same fragment
        for bond in atom.bonds:
            if bond in remaining_bonds or bond.mirror in remaining_bonds:
                assign_to_fragment(bond.a2, frag_id)
    num_frag_id = 0
    for atom in molfile.atoms:
        # We'll skip all but the first atom in each fragment as all atoms
        # in the fragment will be assigned by assign_to_fragment()
        if atom.fragment_id == 0:
            num_frag_id += 1
            assign_to_fragment(atom, num_frag_id)
    #for atom in molfile.atoms:
    #    print atom.name, atom.fragment_id
    # Assert that all atoms have been assigned to a fragment
    assert(len([a for a in molfile.atoms if a.fragment_id == 0]) == 0)
    # More than 9 fragments will break our current residue naming scheme,
    # which uses two letters plus a number from 1 to 9.
    if num_frag_id > 9:
        raise ValueError("More than 9 ligand fragments!")
    # No atom may have more than 1 connection outside the residue.
    # No two fragments may have more than one connection.
    frag_frag_conns = set()
    for atom in molfile.atoms:
        if len(atom.conn_bonds) > 1:
            pass #raise ValueError("Cannot create more than one connection at same atom (%s)" % atom.name)
        # Only check one permutation b/c the other gets caught going backwards along the same bond.
        for conn in atom.conn_bonds:
            pair1 = (conn.a1.fragment_id, conn.a2.fragment_id)
            #pair2 = (conn.a2.fragment_id, conn.a1.fragment_id)
            if pair1 in frag_frag_conns:
                #raise ValueError("Cannot create multiple connections between fragments %i and %i" % pair1)
                #print "WARNING: Multiple connections between fragments (%i and %i) **NOT CURRENTLY SUPPORTED** by Rosetta!" % pair1
                pass # not a concern anymore: is now supported by Mini (I think)
            frag_frag_conns.add(pair1)
            #frag_frag_conns.add(pair2)
    # Fragments should be about the same size as other residues
    for frag_id in range(1,num_frag_id+1):
        frag_atoms = [a for a in molfile.atoms if a.fragment_id == frag_id]
        num_atoms = len(frag_atoms)
        num_heavy_atoms = len([a for a in frag_atoms if not a.is_H])
        num_rot_bonds = len([b for b in molfile.bonds if b.a1.fragment_id == frag_id and b.a2.fragment_id == frag_id and b.can_rotate])
        if num_atoms < 3:
            # Mini-Rosetta atomtree requires at least 3 points to establish a coordinate system.
            print "Fragment %i: %s" % (frag_id, [a.name for a in frag_atoms])
            raise ValueError("Fragment %i has %i atoms; merge with another fragment or add virtual atoms to make 3 total" % (frag_id, num_atoms))
        if not (7 <= num_atoms <= 24):
            print "WARNING: fragment %i has %i total atoms including H; protein residues have 7 - 24 (DNA: 33)" % (frag_id, num_atoms)
        if not (4 <= num_heavy_atoms <= 22):
            print "WARNING: fragment %i has %i non-H atoms; protein residues have 4 - 14 (DNA: 22)" % (frag_id, num_heavy_atoms)
        if num_rot_bonds > 4:
            print "WARNING: fragment %i has %i rotatable bonds; protein residues have 0 - 4" % (frag_id, num_rot_bonds)
    print "Average %.1f atoms (%.1f non-H atoms) per fragment" % (
        float(len(molfile.atoms)) / float(num_frag_id), float(len([a for a in molfile.atoms if not a.is_H])) / float(num_frag_id))
    # Average stats tabulated by IWD from Richardson's Top500 database
    print "(Proteins average 15.5 atoms (7.8 non-H atoms) per residue)"
    return num_frag_id
#}}}
def build_fragment_trees(molfile): #{{{
    '''Assigns a root atom for each fragment and parents and children.'''
    # Assign root atoms based on instructions in the molfile
    for line in molfile.footer:
        # Standard MDL style is with a space, but KWK has used "MROOT" in the past.
        if line.startswith("M ROOT"):  molfile.atoms[ int(line.split()[2]) - 1 ].is_root = True
        elif line.startswith("MROOT"): molfile.atoms[ int(line.split()[1]) - 1 ].is_root = True
        elif line.startswith("M NBR"):  molfile.atoms[ int(line.split()[2]) - 1 ].is_nbr = True
    # Clear out previously assigned parents / children, if any:
    for a in molfile.atoms:
        a.parent = None
        a.children = []
    for frag_id in set([a.fragment_id for a in molfile.atoms]):
        # If we want to have a default way of choosing the root atom, this is the place:
        root_atoms = [a for a in molfile.atoms if a.fragment_id == frag_id and a.is_root]
        if len(root_atoms) == 0:
            print "WARNING:  no root atom specified, using auto-selected NBR atom instead."
            (nbr, nbr_dist) = choose_neighbor_atom(molfile, frag_id)
            nbr.is_root = True
            root_atoms = [nbr]
        elif len(root_atoms) != 1:
            raise ValueError("You must have no more than one 'M ROOT' record in the molfile per ligand fragment")
        root_atom = root_atoms[0]
        # Protein residues appear to go depth first, so that all chi angles ride on each other.
        # Depth first assignment -- leads to very deep trees
        def tree_dfs(parent):
            # Want to visit non-H children first
            tmp_children = [b.a2 for b in parent.bonds]
            tmp_children.sort(lambda a,b: cmp(a.is_H, b.is_H))
            for child in tmp_children:
                if child.fragment_id != parent.fragment_id: continue
                if child.parent is not None or child.is_root: continue
                child.parent = parent
                parent.children.append(child)
                tree_dfs(child)
        tree_dfs(root_atom)
        ## Breadth first assigment -- minimizes depth of the tree
        #bfs_list = [root_atom]
        #while len(bfs_list) > 0:
        #    # pop first item
        #    parent = bfs_list[0]
        #    del bfs_list[0]
        #    for bond in parent.bonds:
        #        child = bond.a2
        #        if child.fragment_id != parent.fragment_id: continue
        #        if child.parent is not None or child.is_root: continue
        #        child.parent = parent
        #        parent.children.append(child)
        #        bfs_list.append(child)
        #    # sort heavy atom children before hydrogens
        #    parent.children.sort(lambda a,b: cmp(a.is_H, b.is_H))
    # Every atom should have a parent OR be a root
    assert(len([a for a in molfile.atoms if not a.is_root and a.parent is None]) == 0)
#}}}
def assign_internal_coords(molfile): #{{{
    '''Sets up stubs/input_stubs and d,theta,phi for all atoms.'''
    for frag_id in set([a.fragment_id for a in molfile.atoms]):
        root_atoms = [a for a in molfile.atoms if a.fragment_id == frag_id and a.is_root]
        assert(len(root_atoms) == 1)
        root_atom = root_atoms[0]
        # Assign stub atoms so we can calculate our internal coords
        # Root atom doesn't use input stubs -- has no internal coords
        def assign_stubs(me):
            # Find stub atoms for me.  Must store them b/c defined recursively.
            # Logic copied from core::kinematics::tree::Atom_, BondedAtom_, JumpAtom_ (*.hh)
            # Only the root atom is a jump in these cases, simplifying things somewhat.
            if me.is_root:
                me.stub1 = me
                me.stub2 = me.children[0] # first child
                if len(me.stub2.children) > 0:  me.stub3 = me.stub2.children[0] # first child of first child
                else:                           me.stub3 = me.children[1] # second child
            else:
                me.stub1 = me
                me.stub2 = me.parent
                me.stub3 = me.parent.stub2
                # Special case for first child of the root
                if me.parent.is_root and me.stub3 == me:
                    me.stub3 = me.parent.stub3
            #print "stubs", [x.name for x in (me, me.stub1, me.stub2, me.stub3)]
            # Assign input stubs to children and calculate internal coords
            prev_sibling = me.stub3
            parent = me # rename to make logic clearer
            # Have to store input_stub atoms b/c they're written to params file
            for child in parent.children:
                child.input_stub1 = parent.stub1
                child.input_stub2 = parent.stub2
                child.input_stub3 = prev_sibling # for first child, this is parent.stub3
                # Special case for second child of the root
                if parent.is_root and prev_sibling == parent.stub2:
                    #print "activate second child case! stub3 =", parent.stub3.name
                    child.input_stub3 = parent.stub3
                #print "input_stubs", [x.name for x in (child, child.parent, child.input_stub1, child.input_stub2, child.input_stub3)]
                # Now actually calculate spherical internal coordinates
                child.d, child.theta, child.phi = calc_internal_coords(child, child.input_stub1, child.input_stub2, child.input_stub3)
                # Recursive update of child's children
                assign_stubs(child)
                # Child is now previous sibling for next child in this loop
                prev_sibling = child
        # end assign_stubs()
        assign_stubs(root_atom)
        # Root has to have dummy input stub atoms to fill space in params file
        root_atom.input_stub1 = root_atom.stub1
        root_atom.input_stub2 = root_atom.stub2
        root_atom.input_stub3 = root_atom.stub3
        # Root has dummy values for d, theta, phi
        root_atom.d     = 0.0
        root_atom.theta = 0.0
        root_atom.phi   = 0.0
    # end loop over all fragments
#}}}
def calc_internal_coords(child, input_stub1, input_stub2, input_stub3): #{{{
    '''Returns (d, theta, phi) for a point given it's three input stub points.'''
    #print "calc_internal_coords", [x.name for x in (child, input_stub1, input_stub2, input_stub3)]
    # Now actually calculate spherical internal coordinates
    # There's some weird "flip_stub" logic for first child of the root
    # iff theta is to be kept fixed (BondedAtom.cc) but I don't think it matters here.
    #
    # parent == parent.stub1 == child.input_stub1, always
    # (except for CONNECT atoms, where we use different stubs!!)
    d = r3.distance(child, input_stub1)
    if d < 1e-2: # very small d;  theta, phi don't make sense
        print "WARNING: very small d=%f for %s" % (d, child.name)
        theta = 0.0
        phi = 0.0
    else:
        theta = r3.angle( r3.from_to(input_stub2,input_stub1), r3.from_to(input_stub1,child) )
        if theta < 1e-2 or theta > 180 - 1e-2:
            # This always happens for first child of root:
            #print "WARNING: nearly parallel theta=%f for %s" % (theta, child.name)
            phi = 0.0
        else:
            phi = r3.dihedral(child, input_stub1, input_stub2, input_stub3)
    return (d, theta, phi)
#}}}
@memoize
def choose_neighbor_atom(molfile, frag_id, nbr_atom=None): #{{{
    '''You may specify the desired nbr_atom manually if you just want to compute the nbr_dist.'''
    atoms = [atom for atom in molfile.atoms if atom.fragment_id == frag_id]
    bonds = [bond for bond in molfile.bonds if bond.a1.fragment_id == frag_id and bond.a2.fragment_id == frag_id]
    # If the entire fragment is rigid, max distance between any two atoms
    # is their Cartesian distance.
    # If not, an upper bound on the max distance between two atoms is the sum
    # of the bond lengths along the path between them that minimizes said sum.
    # We can make this bound tighter by doing Cartesian distances within a
    # rigid substructure and only summing bonds between rigid units.
    # We implement this by making all atoms within a rigid unit
    # pseudo-neighbors of each other, joined together by pseudo-bonds.
    #
    # Create an N x N array of all-against-all distances.
    # Syntax is awkward and Python lacks a reliable +Inf, so we use 1e100.
    na = len(atoms) # Number of Atoms
    all_all_dist = [ [1e100] * na for i in range(na) ]
    num_rot_bonds = len([b for b in bonds if b.can_rotate])
    if num_rot_bonds == 0: # easy case -- rigid fragment
        for i in range(0,na):
            all_all_dist[i][i] = 0
            for j in range(i+1,na):
                d = r3.distance( atoms[i], atoms[j] )
                all_all_dist[i][j] = d
                all_all_dist[j][i] = d
    else: # hard case -- flexible fragment
        # This will be queried several times so better to precompute it:
        nbrs = dict([ (a,set()) for a in atoms ])
        for a in atoms:
            # a's neighbors:  bonded and in same fragment
            nbrs[a].update([b.a2 for b in a.bonds if b.a2.fragment_id == frag_id])
            # a's pseudo-neighbors:  in same rigid unit and same fragment
            nbrs[a].update([a2 for a2 in atoms if a2.rigid_id == a.rigid_id])
        # The Floyd-Warshall version is significantly faster (~2x)
        #for i in range(0,na):
        #    all_all_dist[i] = dijkstra(
        #        start = atoms[i],
        #        nodes = atoms,
        #        nbr = lambda a: nbrs[a],
        #        dist = r3.distance
        #    )
        all_all_dist = floyd_warshall(
                nodes = atoms,
                nbr = lambda a: nbrs[a],
                dist = r3.distance
            )
    # The "best" neighbor atom is the one with
    # the smallest {max distance to all other atoms}.
    max_dist = [ max(row) for row in all_all_dist ]
    # If not specified as an argument, did the user request a specific neighbor atom?
    if nbr_atom is None:
        nbr_atoms = [a for a in atoms if a.is_nbr]
        assert len(nbr_atoms) <= 1
        if len(nbr_atoms) > 0: nbr_atom = nbr_atoms[0]
    if nbr_atom is None:
        # For a robust atom tree, however (especially converting centroid -> fullatom),
        # we want a heavy atom with 2+ heavy atom neighbors:
        one_is_ok = False
        for i, atom in enumerate(atoms):
            num_heavy_bonds = len([b for b in atom.heavy_bonds if b.a2.fragment_id == frag_id])
            if atom.is_H or num_heavy_bonds < 2: max_dist[i] = 1e100
            else: one_is_ok = True
        if not one_is_ok: raise ValueError("No acceptable neighbor atoms in molecule!")
        best_idx = argmin(max_dist)
    else: best_idx = atoms.index(nbr_atom)
    return (atoms[best_idx], max_dist[best_idx])
#}}}
def floyd_warshall(nodes, nbr, dist): #{{{
    '''Computes the shortest path from each node to all others.
    nodes - the set of all nodes
    nbr() - returns a list of a node's neighbors, nbr(n) = [n's neighbors]
    dist() - returns distance between two nodes that are neighbors, dist(n,m)
    Returns 2D array of shortest distances in same order as input nodes.
    '''
    # Create an N x N array of all-against-all distances.
    # Syntax is awkward and Python lacks a reliable +Inf, so we use 1e100.
    N = len(nodes) # Number of Atoms
    range_N = range(N)
    d = [ [1e100] * N for i in range_N ]
    # Initialization: distance between directly connected nodes
    for j in range_N:
        nbrs_j = nbr(nodes[j])
        for i in range_N:
            if nodes[i] in nbrs_j:
                d[i][j] = dist(nodes[i], nodes[j])
    # Propagation: expand set of usable intermediate vertices one at a time
    for k in range_N:
        for j in range_N:
            for i in range_N:
                d[i][j] = min(d[i][j], d[i][k]+d[k][j])
    return d
#}}}
def dijkstra(start, nodes, nbr, dist): #{{{
    '''Computes the shortest path from start node to all others.
    start - the node to start from
    nodes - the set of all nodes (includes start)
    nbr() - returns a list of a node's neighbors, nbr(n) = [n's neighbors]
    dist() - returns distance between two nodes that are neighbors, dist(n,m)
    Returns array of shortest distances to start in same order as input nodes.
    '''
    # Not as efficient as it could be, I think, but best I could manage.
    # Tmp objects [dist_from_start,node] sort properly
    DIST = 0; NODE = 1
    queue = [ [1e100,node] for node in nodes ] # 1e100  ~  +Inf
    # Allows lookup of best distance by name
    shortest = dict([ (q[NODE],q) for q in queue ])
    shortest[start][DIST] = 0 # start
    while len(queue) > 0:
        curr_idx = argmin(queue) # on first pass this is start
        curr = queue.pop(curr_idx)[NODE]
        curr_shortest = shortest[curr][DIST]
        for n in nbr(curr):
            new_dist = curr_shortest + dist(curr,n)
            if new_dist < shortest[n][DIST]:
                shortest[n][DIST] = new_dist
    return [ shortest[node][DIST] for node in nodes ]
#}}}
def write_ligand_kinemage(f, molfile): #{{{
    if not isinstance(f, file): f = gz_open(f, 'w')
    f.write("@text\n")
    f.write("View this file with KiNG or Mage from http://kinemage.biochem.duke.edu\n")
    f.write("@kinemage 1\n")
    f.write("@title {%s}\n" % os.path.basename(f.name))
    f.write("@onewidth\n")
    # Element markers
    elem_colors = { "C": "green", "N": "sky", "O": "red", "H": "gray", "S": "yellow", "P": "peach", "F": "bluetint", "CL": "cyan", "BR": "sea", "I": "lilac" }
    f.write("@balllist {element balls} color= gray radius= 0.1\n")
    for a in molfile.atoms: f.write("{%s} %s %.3f %.3f %.3f\n" % (a.elem, elem_colors.get(a.elem.upper(), "hotpink"), a.x, a.y, a.z))
    # All bonds, including which can rotate
    f.write("@vectorlist {all bonds} color= gray width= 1\n")
    for b in molfile.bonds:
        f.write("{%s}P %.3f %.3f %.3f\n" % (b.a1.name, b.a1.x, b.a1.y, b.a1.z))
        f.write("{%s}L %.3f %.3f %.3f\n" % (b.a2.name, b.a2.x, b.a2.y, b.a2.z))
    f.write("@vectorlist {rotatable bonds} color= white width= 4\n")
    for b in molfile.bonds:
        if not b.can_rotate: continue
        f.write("{%s}P %.3f %.3f %.3f\n" % (b.a1.name, b.a1.x, b.a1.y, b.a1.z))
        f.write("{%s}L %.3f %.3f %.3f\n" % (b.a2.name, b.a2.x, b.a2.y, b.a2.z))
    # Atom labels
    ai = index_atoms(molfile.atoms) # Atom order in the molfile
    f.write("@labellist {atom indices} color= white off\n")
    for a in molfile.atoms: f.write("{%i} %.3f %.3f %.3f\n" % (ai[a], a.x, a.y, a.z))
    f.write("@labellist {original names} color= white off\n")
    for a in molfile.atoms: f.write("{%s} %.3f %.3f %.3f\n" % (a.orig_name, a.x, a.y, a.z))
    f.write("@labellist {atom names} color= white off\n")
    for a in molfile.atoms: f.write("{%s} %.3f %.3f %.3f\n" % (a.name, a.x, a.y, a.z))
    f.write("@labellist {Rosetta types} color= white\n")
    for a in molfile.atoms: f.write("{%s} %.3f %.3f %.3f\n" % (a.ros_type, a.x, a.y, a.z))
    f.write("@labellist {MM types} color= white off\n")
    for a in molfile.atoms: f.write("{%s} %.3f %.3f %.3f\n" % (a.mm_type, a.x, a.y, a.z))
    f.write("@labellist {partial charges} color= white off\n")
    for a in molfile.atoms: f.write("{%.2f} %.3f %.3f %.3f\n" % (a.partial_charge, a.x, a.y, a.z))
    # Each fragment and its atom tree, distinguished by color
    colors = ['deadwhite', 'purple', 'blue', 'sky', 'cyan', 'sea', 'green', 'lime', 'yellow', 'gold', 'orange', 'red']
    frag_ids = list(set([a.fragment_id for a in molfile.atoms]))
    frag_ids.sort()
    for frag_id in frag_ids:
        color = colors[ frag_id % len(colors) ]
        f.write("@group {frag %i} master= {fragments}\n" % frag_id)
        # Find and mark root atom for this fragment (should only be one)
        f.write("@balllist {root atom} color= %s radius= 0.2 master={root atoms}\n" % color)
        root_atoms = [a for a in molfile.atoms if a.fragment_id == frag_id and a.is_root]
        for a in root_atoms: f.write("{%s frag %i} %.3f %.3f %.3f\n" % (a.name, a.fragment_id, a.x, a.y, a.z))
        (nbr, nbr_dist) = choose_neighbor_atom(molfile, frag_id)
        f.write("@ringlist {nbr atom} color= %s radius= %.3f master={nbr atoms} off\n" % (color, nbr_dist))
        f.write("{%s frag %i} %.3f %.3f %.3f\n" % (nbr.name, nbr.fragment_id, nbr.x, nbr.y, nbr.z))
        f.write("{%s frag %i} r=0.3 %.3f %.3f %.3f\n" % (nbr.name, nbr.fragment_id, nbr.x, nbr.y, nbr.z))
        f.write("@arrowlist {atom tree} color= %s master= {atom trees}\n" % color) # can specify width, angle, radius
        frag_atoms = [a for a in molfile.atoms if a.fragment_id == frag_id]
        for a in frag_atoms:
            if a.parent is None: continue
            a1 = a.parent; a2 = a
            f.write("{%s}P %.3f %.3f %.3f\n" % (a1.name, a1.x, a1.y, a1.z))
            f.write("{%s}L %.3f %.3f %.3f\n" % (a2.name, a2.x, a2.y, a2.z))
        # Connection arrows are half length and meet in middle of bond
        f.write("@arrowlist {connections} color= %s master= {connections}\n" % color) # can specify width, angle, radius
        for a1 in frag_atoms:
            for b in a1.conn_bonds:
                #a2 = r3.midpoint(b.a1, b.a2)
                # so they don't *quite* meet in the middle (45% across):
                a2 = r3.add( r3.mult(b.a1,0.55), r3.mult(b.a2,0.45) )
                f.write("{%s}P %.3f %.3f %.3f\n" % (a1.name, a1.x, a1.y, a1.z))
                f.write("{%s}L %.3f %.3f %.3f\n" % ("midpt", a2.x, a2.y, a2.z))
    f.close()
#}}}
def write_param_file(f, molfile, name, frag_id, base_confs, max_confs, amino_acid=None,long_names=False,write_pdb_rotamers=False): #{{{
    '''Writes a Molfile object to a file.
    f may be a file name or file handle.
    base_confs is the number of heavy-atom conformations generated by e.g. Omega
    max_confs is the maximum number of conformations desired after adding proton rotation
        The total number of confs may still be larger than max_confs,
        but at least we'll skip -ex# extra sampling of proton chis.
    '''
    close_file = False
    if not isinstance(f, file):
        f = gz_open(f, 'w')
        close_file = True
    full_name = name # keep the untruncated name so we can use it with the long_names option
    if frag_id == 1 and len(name) > 2: name = "%3.3s" % name
    else: name = "%2.2s%1i" % (name, frag_id)
    name_string = "NAME %s\n"
    if long_names:
        f.write(name_string % full_name)
    else:
        f.write(name_string % name)
    if amino_acid is None:
        f.write("IO_STRING %3.3s Z\n" % name) # 'Z' is the agreed-upon one letter code for ligands in Mini
        f.write("TYPE LIGAND\n")
        f.write("AA UNK\n")
    else:
        f.write("IO_STRING %3.3s X\n" % name)
        #f.write("AA %3.3s\n" % amino_acid)
        f.write("AA UNK\n") # oops, looks like you want aa_unk even for modified amino acids
        f.write("TYPE POLYMER\n")
        f.write("LOWER_CONNECT N\n")
        f.write("UPPER_CONNECT C\n")
        f.write("PROPERTIES PROTEIN\n")
        f.write("FIRST_SIDECHAIN_ATOM CB\n")
        f.write("ACT_COORD_ATOMS CB END\n")
    atoms = [atom for atom in molfile.atoms if atom.fragment_id == frag_id]
    bonds = [bond for bond in molfile.bonds if bond.a1.fragment_id == frag_id and bond.a2.fragment_id == frag_id]
    root_atoms = [atom for atom in atoms if atom.is_root]
    assert(len(root_atoms) == 1)
    root_atom = root_atoms[0]
    # ORDER that atoms appear seems to imply the order of the atom tree,
    # rather than the order of the ICOOR_INTERNAL records (?)
    #for atom in atoms:
    #    f.write("ATOM %-4s %-4s %-4s %.2f\n" % (atom.name, atom.ros_type, atom.mm_type, atom.partial_charge))
    # So instead we write atoms depth-first, starting from root.
    def write_atoms(atom):
        f.write("ATOM %-4s %-4s %-4s %.2f\n" % (atom.name, atom.ros_type, atom.mm_type, atom.partial_charge))
        for a2 in atom.children: write_atoms(a2)
    write_atoms(root_atom)
    for bond in bonds:
        f.write("BOND_TYPE %-4s %-4s %-4s\n" % (bond.a1.name, bond.a2.name, bond.order))
    # Define chi angles
    # Iterating over the bonds is non-trivial and we need multiple passes, so we define a generator:
    def rot_bond_iter(bonds):
        '''Yields the tuples (bond, a, b, c, d) for each rotatable bond and the four chi atoms.'''
        for bond in bonds:
            if not bond.can_rotate: continue
            # define atoms a-b-c-d with d closest to tree root
            if bond.a1 == bond.a2.parent:
                b = bond.a2
                c = bond.a1
            elif bond.a2 == bond.a1.parent:
                b = bond.a1
                c = bond.a2
            else: raise ValueError("Rotatable bond %s - %s not in atom tree!" % (bond.a1.name, bond.a2.name))
            # Child MUST have a child for this to have counted as rotatable.
            # Might as well use the first one, it's more likely to be a heavy atom.
            # ... unless child is at end of tree and has no other children except
            #                                   |    |                  |
            # the CONN atom, e.g. an ether:  -> D -> C -> B(oxygen) -|> A(next frag)
            #                                   |    |                  |
            # This seems likely to blow Rosetta's mind, so let's not allow it:
            if len(b.children) < 1:
                raise ValueError("Can't define chi angle: no children for %s.  Don't split ether bonds, OK?  Try --no-param for debugging." % b.name)
            a = b.children[0]
            if c.parent is not None:
                d = c.parent
            else:
                # If c is root atom, have to choose a different anchor
                # c MUST have 2+ children for this bond to have been called rotatable
                # (if it doesn't, we made a logical error somewhere...)
                # So, use first child of c that's not b:
                # ... again, unless C is root and the O of an ether bond (no other kids)
                d = [k for k in c.children if k != b][0]
            yield (bond, a, b, c, d)
    def is_sp2_proton(a, b, c, d):
        '''Crude guestimate of H-bonding geometry'''
        # i.e., assume sp2 if the H's grandparent has anything other than single bonds to it.
        return ((a.is_H and len([bnd for bnd in c.heavy_bonds if bnd.order != Bond.SINGLE]) > 0)
             or (d.is_H and len([bnd for bnd in b.heavy_bonds if bnd.order != Bond.SINGLE]) > 0))
    # Do proton chi's first so we can use -ex1, -ex2, etc
    sorted_bonds = list(bonds) # make a copy then sort in place
    sorted_bonds.sort(key=lambda b: b.is_proton_chi, reverse=True)
    num_H_confs = base_confs
    for bond, a, b, c, d in rot_bond_iter(sorted_bonds):
        if bond.is_proton_chi:
            if is_sp2_proton(a, b, c, d): num_H_confs *= 6
            else: num_H_confs *= 9
    if num_H_confs > max_confs:
        print "WARNING: skipping extra samples for proton chis; would give %i conformers" % num_H_confs
    num_chis = 0
    for bond, a, b, c, d in rot_bond_iter(sorted_bonds):
        num_chis += 1
        # "Direction" of chi definition appears not to matter,
        # but we follow the convention of amino acids (root to tip)
        f.write("CHI %i %-4s %-4s %-4s %-4s\n" % (num_chis, d.name, c.name, b.name, a.name))
        if bond.is_proton_chi:
            # Only expand proton chis with extra samples if doing so won't generate tens of thousands of confs.
            if num_H_confs <= max_confs: extra = "1 20"
            else: extra = "0"
            if is_sp2_proton(a, b, c, d):
                f.write("PROTON_CHI %i SAMPLES 2 0 180 EXTRA %s\n" % (num_chis, extra))
            else:
                f.write("PROTON_CHI %i SAMPLES 3 60 -60 180 EXTRA %s\n" % (num_chis, extra))
    # Assign numbers to the connection points in this fragment
    num_conns = 0
    conn_nums = {}
    for a in atoms:
        for b in a.conn_bonds:
            num_conns += 1
            conn_nums[b] = num_conns
            if b.can_rotate: rot_flag = "CAN_ROTATE"
            else: rot_flag = "NO_ROTATE"
            f.write("CONNECT %-4s %s #CONN%i\n" % (a.name, rot_flag, num_conns));
    if amino_acid is None:
        (nbr, nbr_dist) = choose_neighbor_atom(molfile, frag_id)
    else:
        cb_atoms = [a for a in atoms if a.name == " CB "]
        assert len(cb_atoms) == 1
        (nbr, nbr_dist) = choose_neighbor_atom(molfile, frag_id, cb_atoms[0])
        assert nbr.name == " CB "
    f.write("NBR_ATOM %s\n" % nbr.name)
    f.write("NBR_RADIUS %f\n" % nbr_dist)
    # Formal charges
    for a in atoms:
        if a.formal_charge:
            f.write("CHARGE %s FORMAL %d\n" % (a.name, a.formal_charge))
    # Convention seems to be a depth-first traversal from the root.
    # I don't know whether this matters, but it's the easy way anyhow.
    def write_icoords(a):
        f.write("ICOOR_INTERNAL   %-4s %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n"
            % (a.name, a.phi, a.theta, a.d, a.input_stub1.name, a.input_stub2.name, a.input_stub3.name));
        # Have to re-create some of the internal coord logic here for connection points
        # The bad thing about the current arrangement is that CONN points
        # may be come after hydrogens and use them as reference points...
        prev_sibling = a.stub3
        for child in a.children:
            write_icoords(child)
            prev_sibling = child
        for conn in a.conn_bonds:
            inp_stub1 = a.stub1
            inp_stub2 = a.stub2
            inp_stub3 = prev_sibling
            # Special case when connection is 2nd child of root
            if inp_stub3 == inp_stub2:
                assert(a.is_root)
                inp_stub3 = a.stub3
            d, theta, phi = calc_internal_coords(conn.a2, inp_stub1, inp_stub2, inp_stub3)
            f.write("ICOOR_INTERNAL  CONN%1i %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n"
                % (conn_nums[conn], phi, theta, d, inp_stub1.name, inp_stub2.name, inp_stub3.name));
            # This breaks.  Likewise, using the previous CONNx point breaks. Empirically,
            # the right answer seems to be to keep using the last child atom as inp_stub3.
            #prev_sibling = conn.a2 # if we allowed multiple connects from one parent
    write_icoords(root_atom)
    if amino_acid is not None:
        f.write("# These lines stolen from ALA.params\n")
        f.write("# I *believe* they will override any duplicate definitions above\n")
        f.write("ICOOR_INTERNAL  UPPER  149.999985   63.800007    1.328685   C     CA    N  \n")
        f.write("ICOOR_INTERNAL    O   -180.000000   59.200005    1.231015   C     CA  UPPER\n")
        f.write("ICOOR_INTERNAL  LOWER -150.000000   58.300003    1.328685   N     CA    C  \n")
        f.write("ICOOR_INTERNAL    H   -180.000000   60.849998    1.010000   N     CA  LOWER\n")
    if write_pdb_rotamers and os.path.exists("%s_conformers.pdb" % full_name):
        f.write("PDB_ROTAMERS %s_conformers.pdb\n" % full_name)
    #
    # XXX-FIXME:  still might need PROPERTIES, FIRST_SIDECHAIN_ATOM, ACT_COORD_ATOMS
    #
    if close_file: f.close()
#}}}
def write_ligand_pdb(f, molfile_tmpl, molfile_xyz, resname, ctr=None, chain_id='X'): #{{{
    '''Writes a PDB file with a series of residues representing one ligand conformation.
    The topology (atom names, fragment divisions, etc) are taken from molfile_tmpl,
    while the actual XYZ coordinates are taken from molfile_xyz.
    resname provides the first two characters of the residue name.
    f may be a file name or file handle.'''
    if not isinstance(f, file): f = gz_open(f, 'w')
    # If ctr is set, make it an offset vector for recentering ligands
    if ctr is not None:
        curr_ctr = r3.centroid([a for a in molfile_xyz.atoms if not a.is_H])
        ctr = r3.sub(ctr, curr_ctr)
    else: ctr = r3.Triple(0,0,0)
    atom_num = 0
    frag_ids = list(set([a.fragment_id for a in molfile_tmpl.atoms]))
    frag_ids.sort()
    for frag_id in frag_ids:
        ai = index_atoms(molfile_tmpl.atoms) # 1-based index
        atoms = [a for a in molfile_tmpl.atoms if a.fragment_id == frag_id]
        
        #loop through atoms and make sure the atoms are in the right order before outputting anything
        #this will avoid the possibility of writing half a conformer and then dying
        for atom_tmpl in atoms:
            atom_xyz = molfile_xyz.atoms[ ai[atom_tmpl]-1 ]
            # Somewhere between version 2.2.1 and 2.3.1, Omega started writing the atoms for
            # the -includeInput conformer in a different order from all the other, generated conformers.
            # (According to OpenEye, this specifically affects PDB inputs with MOL2 outputs.)
            # This totally screws the output of this script, as all the PDB rotamers have atoms in the wrong places.
            # This check will catch the problem sometimes, though obviously not always:
            if atom_tmpl.elem != atom_xyz.elem and not atom_tmpl.is_virtual:
                f.write("TER"+(" "*77)+"\n")
                f.close()
                raise ValueError("Atom order mismatch between first and subsequent conformer: %i    %s    %s" % (ai[atom_tmpl], atom_tmpl, atom_xyz))
        
        for atom_tmpl in atoms:
            atom_xyz = molfile_xyz.atoms[ ai[atom_tmpl]-1 ]
            xyz = r3.add(atom_xyz, ctr)
            atom_num += 1
            if len(frag_ids) == 1 and len(resname) > 2:
                f.write("HETATM%5i %-4.4s %3.3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  \n"
                    % (atom_num, atom_tmpl.name, resname,          chain_id, frag_id, xyz.x, xyz.y, xyz.z, 1.0, 20.0, atom_tmpl.elem))
            else:
                f.write("HETATM%5i %-4.4s %2.2s%1i %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  \n"
                    % (atom_num, atom_tmpl.name, resname, frag_id, chain_id, frag_id, xyz.x, xyz.y, xyz.z, 1.0, 20.0, atom_tmpl.elem))
    f.write("TER"+(" "*77)+"\n")
    f.close()
#}}}
def write_fragment_mol2(f, molfile, frag_id): #{{{
    '''Writes a .mol2 file consisting only of atoms from a particular fragment.'''
    m = copy.copy(molfile)
    # Write out only non-hydrogen atoms so that Omega can fill in H on dummies as needed!
    #m.atoms = [ a for a in m.atoms if a.fragment_id == frag_id and not a.is_H ]
    #m.bonds = [ b for b in m.bonds if b.a1.fragment_id == frag_id and b.a2.fragment_id == frag_id and not b.a1.is_H and not b.a2.is_H ]
    m.atoms = [ a for a in m.atoms if a.fragment_id == frag_id ]
    m.bonds = [ b for b in m.bonds if b.a1.fragment_id == frag_id and b.a2.fragment_id == frag_id ]
    write_tripos_mol2(f, [m])
#}}}
def write_all_files(m, molfiles, num_frags, options, suffix=""): #{{{
    '''Returns 0 for success, non-zero error code for fatal error'''
    
    if not options.no_pdb:
        if options.conformers_in_one_file:
            pdb_file = "%s%s.pdb" % (options.pdb,suffix)
            conformer_file = "%s%s_conformers.pdb" % (options.pdb,suffix)
            for i,molfile in enumerate(molfiles):
                if i == 0:
                    file_to_write = pdb_file
                else:
                    file_to_write = open(conformer_file,'a')
                try:
                    write_ligand_pdb(file_to_write, m, molfile, options.name, options.center)
                except ValueError as e:
                    sys.exit(e)
        else:
            for i, molfile in enumerate(molfiles):
                pdb_file = "%s_%04i%s.pdb" % (options.pdb, i+1, suffix)
                if not options.clobber and os.path.exists(pdb_file):
                    print "File %s already exists -- aborting!" % pdb_file
                    print "Use --clobber to overwrite existing files."
                    return 4
                else:
                    # m is used for names, molfile is used for XYZ
                    try:
                        write_ligand_pdb(pdb_file, m, molfile, options.name, options.center)
                    except ValueError as e:
                        if options.skip_bad_conformers:
                            print "Skipping Bad Conformers:",e
                            os.remove(pdb_file)
                        else:
                            sys.exit(e)   
                    print "Wrote PDB file %s" % pdb_file
    if not options.no_param:
        for i in range(num_frags):
            if num_frags == 1: param_file = "%s%s.params" % (options.pdb, suffix)
            else: param_file = "%s%i%s.params" % (options.pdb, i+1, suffix)
            if not options.clobber and os.path.exists(param_file):
                print "File %s already exists -- aborting!" % param_file
                print "Use --clobber to overwrite existing files."
                return 2
            else:
                if len(molfiles) > 1:
                    write_param_file(param_file, m, options.name, i+1, len(molfiles), options.max_confs, options.amino_acid, long_names=options.long_names, write_pdb_rotamers=options.conformers_in_one_file)
                else:
                    write_param_file(param_file, m, options.name, i+1, len(molfiles), options.max_confs, options.amino_acid, long_names=options.long_names, write_pdb_rotamers=False)
                print "Wrote params file %s" % param_file
    if options.kinemage is not None:
        kin_file = options.kinemage
        dot = kin_file.rfind(".")
        if dot != -1: kin_file = kin_file[:dot] + suffix + kin_file[dot:]
        if not options.clobber and os.path.exists(kin_file):
            print "File %s already exists -- aborting!" % kin_file
            print "Use --clobber to overwrite existing files."
            return 3
        else:
            write_ligand_kinemage(kin_file, m)
            print "Wrote kinemage file %s" % kin_file

    if num_frags > 1:
        for i in range(num_frags):
            mol2_file = "%s%i%s.mol2" % (options.pdb, i+1, suffix)
            if not options.clobber and os.path.exists(mol2_file):
                print "File %s already exists -- aborting!" % mol2_file
                print "Use --clobber to overwrite existing files."
                return 5
            else:
                write_fragment_mol2(mol2_file, m, i+1)
                print "Wrote mol2 file %s" % mol2_file
    
    if options.extra_torsion_output:
        tor_file = "%s%s.tors" % (options.pdb, suffix)
        f = open(tor_file,'w')
        def expand_atoms_by_one_bond(bonds, atoms):
            atoms_current = [a for a in atoms]
            for b in bonds:
                if b.a1 in atoms_current or b.a2 in atoms_current:
                    if b.a1 not in atoms: atoms.append(b.a1)
                    if b.a2 not in atoms: atoms.append(b.a2)
        for b in m.bonds:
            #if (b.order == Bond.AROMATIC):
            if (b.order == Bond.DOUBLE and b.a1.elem in ["C"]) or (b.order == Bond.AROMATIC):
                conjugated_atoms = [b.a1, b.a2]
                expand_atoms_by_one_bond(m.bonds,conjugated_atoms)
                #expand_atoms_by_one_bond(m.bonds,conjugated_atoms)

                for i in range(2,len(conjugated_atoms)):
                    for j in range(i+1,len(conjugated_atoms)):
                        f.write("%s %s %s %s %s  0  80  2\n"%(options.name, conjugated_atoms[i].name, b.a1.name, b.a2.name, conjugated_atoms[j].name))
        f.close()
        
    return 0
#}}}
def main(argv): #{{{
    """
Converts a small molecule in an MDL Molfile with "M SPLT" and "M ROOT"
records into a series of .params residue definition files for Rosetta.
Also writes out the ligand conformation as PDB HETATMs.

If an SD file is given as input instead, the first entry is used for
generating topology / parameter files, and they all are used for
generating PDB-style coordinates in separate, numbered files.
These multiple files can optionally be concatenated into a single file,
which can then be specified with an additional PDB_ROTAMERS line in the
.params file to include the extra conformations as ligand rotamers.
Multiple models may also be supplied in MOL2 format, which does not support
M ROOT and M SPLT records but does allow for partial charges.
File type is deduced from the extension.

To divide a ligand into fragments by "breaking" bonds (optional):
M SPLT atom_no1 atom_no2

To specify a neighbor atom for a ligand fragment (optional):
M NBR atom_no

To specify a root atom for a ligand fragment (optional):
M ROOT atom_no

The "M" records (M SPLT, M NBR, M ROOT) can alternatively be specified in
a separate control file, which can be used with MOL2 format files.

Note that for ligands with multiple rotamers, Rosetta overlays the ligands
based on the neighbor atom (not the root atom), such that the position of the
neighbor atom and the orientation of the atoms bonded to the neighbor atom is
the same. When using ligand rotamers, it is recommended to confirm that the
neighbor atom falls in an appropriate position.

Expects that the input ligand has already had aromaticity "perceived",
i.e. that it contains aromatic bonds rather than alternating single and double
bonds (Kekule structure).

Optionally writes a kinemage graphics visualization of the atom tree,
neighbor atom selection, fragments, etc -- very helpful for debugging
and for visualizing exactly what was done to the ligand.
    """ # Preformatted
    parser = OptionParser(usage="usage: %prog [flags] { INPUT.mol | INPUT.sdf | INPUT.mol2 }", formatter=PreformattedDescFormatter())
    parser.set_description(main.__doc__)
    # parser.add_option("-short", ["--long"],
    #   action="store|store_true|store_false",
    #   default=True|False|...
    #   type="string|int|float",
    #   dest="opt_name",
    #   help="store value in PLACE",
    #   metavar="PLACE",
    # )
    parser.add_option("-n", "--name",
        default="LG",
        help="name ligand residues NM1,NM2,... instead of LG1,LG2,...",
        metavar="NM"
    )
    parser.add_option("-p", "--pdb",
        default=None, # same as --name, see below
        help="prefix for PDB file names",
        metavar="FILE"
    )
    parser.add_option("-c", "--centroid",
        default=False,
        action="store_true",
        help="write files for Rosetta centroid mode too"
    )
    parser.add_option("--center",
        default=None, # same as --name, see below
        help="translate output PDB coords to have given heavy-atom centroid",
        metavar="X,Y,Z"
    )
    parser.add_option("-m", "--max-confs",
        default=5000, # 400 (default Omega max) * 9 (one sp3 H with -ex1) = 3600
        type="int",
        help="don't expand proton chis if above this many total confs",
        metavar="MAX"
    )
    parser.add_option("-k", "--kinemage",
        default=None,
        help="write ligand topology to FILE",
        metavar="FILE"
    )
    parser.add_option("-a", "--amino-acid",
        default=None,
        help="set up params file for modified amino acid; .mol2 only; edit chis afterward.  Implies --keep-names.",
        metavar="ALA"
    )
    parser.add_option("--clobber",
        default=False,
        action="store_true",
        help="overwrite existing files"
    )
    parser.add_option("--no-param",
        default=False,
        action="store_true",
        help="skip writing .params files (for debugging)"
    )
    parser.add_option("--no-pdb",
        default=False,
        action="store_true",
        help="skip writing .pdb files (for debugging)"
    )
    parser.add_option("--extra_torsion_output",
        default=False,
        action="store_true",
        help="writing additional torsion files"
    )
    parser.add_option("--keep-names",
        default=False,
        action="store_true",
        help="leaves atom names untouched except for duplications"
    )
    parser.add_option("--long-names",
        default=False,
        action="store_true",
        help="if specified name is longer than 3 letters, keep entire name in param NAME field instead of truncating"
    )
    parser.add_option("--recharge",
        type="int",
        help="ignore existing partial charges, setting total charge to CHG",
        metavar="CHG",
    )
    parser.add_option("--m-ctrl",
        default=None,
        help="read additional M control lines from FILE",
        metavar="FILE"
    )
    parser.add_option("--mm-as-virt",
        default=False,
        dest="mm_as_virt",
        help="assign mm atom types as VIRT, rather than X",
        action="store_true"
    )
    parser.add_option("--skip-bad-conformers",
        default=False,
        dest="skip_bad_conformers",
        help="If a conformer has atoms in the wrong order, skip it and continue rather than dying",
        action="store_true"
    )
    parser.add_option("--conformers-in-one-file",
        default=False,
        help="Output 1st conformer to NAME.pdb and all others to NAME_conformers.pdb",
        action="store_true"
    )
    
    
    (options, args) = parser.parse_args(args=argv)
    if options.pdb is None: options.pdb = options.name
    if options.amino_acid is not None: options.keep_names = True

    if len(args) < 1:
        parser.print_help()
        print "Must specify input .mol file!"
        return 1
    elif len(args) == 1:
        infile = args[0]
    else:
        parser.print_help()
        print "Too many arguments!"
        return 1

    ctr = None
    if options.center:
        f = options.center.split(",")
        if len(f) != 3:
            f = options.center.split()
            if len(f) != 3:
                print "Must say -center 'X,Y,Z'"
                return 5
        ctr = r3.Triple( float(f[0]), float(f[1]), float(f[2]) )
        #print "Centering ligands at %s" % ctr

    # There's a very strong order dependence to these function calls:
    # many depend on atom/bond variables set by earlier calls.
    infile_lc = infile.lower()
    if infile_lc.endswith(".mol2") or infile_lc.endswith(".mol2.gz"):
        molfiles = list(read_tripos_mol2(infile, do_find_rings=False))
    elif infile_lc.endswith(".mol") or infile_lc.endswith(".mdl") or infile_lc.endswith(".sdf") or \
			infile_lc.endswith(".mol.gz") or infile_lc.endswith(".mdl.gz") or infile_lc.endswith(".sdf.gz"):
        molfiles = list(read_mdl_sdf(infile, do_find_rings=False))
    else:
        print "Unrecognized file type, must be .mol/.sdf or .mol2!"
        return 6
    # Add additional M_____ control records, if any specified.
    if options.m_ctrl is not None:
        m_ctrl = gz_open(options.m_ctrl,'r')
        try:
            footer = m_ctrl.readlines()
        finally:
            m_ctrl.close()
        for molfile in molfiles:
            molfile.footer.extend(footer)
    m = molfiles[0]
    # Only doing ring detection for the first entry (the only one that needs it)
    # saves a LOT of compute time (~90%) for large input files.
    find_rings(m.bonds)
    # If -center not given, default is to center like first entry
    if ctr is None:
        ctr = r3.centroid([a for a in m.atoms if not a.is_H])
    print "Centering ligands at %s" % ctr
    options.center = ctr

    mark_fragments(m)
    add_fields_to_atoms(m.atoms)
    add_fields_to_bonds(m.bonds)
    find_virtual_atoms(m.atoms)
    if uniquify_atom_names(m.atoms, force=(not options.keep_names)):
        print "Atom names contain duplications -- renaming all atoms."
    for atom in m.atoms: atom.name = pdb_pad_atom_name(atom) # for output convenience
    check_bond_count(m.atoms)
    check_aromaticity(m.bonds)
    check_hydrogens(m.atoms)
    assign_rosetta_types(m.atoms)
    assign_mm_types(m.atoms,options.mm_as_virt)
    assign_centroid_types(m.atoms)
    if options.amino_acid is not None: setup_amino_acid(m.atoms, molfiles)
    net_charge = 0.0
    if options.recharge is not None:
        net_charge = float(options.recharge)
    else:
        #Handle official MDL format formal charges.
        for a in m.atoms:
            net_charge += a.formal_charge
        # KWK's convention -- difference is one space versus two.
        for line in m.footer:
            if line.startswith("M CHG"): net_charge = float(line[5:])
    assign_partial_charges(m.atoms, net_charge, options.recharge is not None)
    assign_rotatable_bonds(m.bonds)
    assign_rigid_ids(m.atoms)
    num_frags = fragment_ligand(m)
    build_fragment_trees(m)
    assign_internal_coords(m)
    #uniquify_atom_names(m.atoms)

    # All-atom output
    if options.centroid:
        ret = write_all_files(m, molfiles, num_frags, options, suffix=".fa")
    else:
        ret = write_all_files(m, molfiles, num_frags, options, suffix="")
    if ret != 0: return ret

    # Centroid mode now
    if options.centroid:
        strip_H(m, lambda a: a.cen_type is None)
        for a in m.atoms: a.ros_type, a.cen_type = a.cen_type, a.ros_type
        build_fragment_trees(m)
        assign_internal_coords(m)
        ret = write_all_files(m, molfiles, num_frags, options, suffix=".cen")
        if ret != 0: return ret

    return 0

if __name__ == "__main__":
    #import cProfile, os
    #i = 1
    #while os.path.exists("profile_%i" % i): i += 1
    #cProfile.run('sys.exit(main(sys.argv[1:]))', "profile_%i" % i)
    sys.exit(main(sys.argv[1:]))

    # Vestigal code for validating automatic atom typing:
    #m = read_mdl_molfile(sys.argv[1])
    #add_fields_to_atoms(m.atoms)
    #assign_rosetta_types(m.atoms)
    #for i, a in enumerate(m.atoms):
    #    err_flag = ""
    #    if a.name.strip() != a.ros_type.strip():
    #        if a.name == "CH1" and a.ros_type.startswith("CH"): pass
    #        else: err_flag = "********" #raise ValueError("typing mismatch!")
    #    print "%3i %4s --> %4s %s" % (i+1, a.name, a.ros_type, err_flag)
#}}}
