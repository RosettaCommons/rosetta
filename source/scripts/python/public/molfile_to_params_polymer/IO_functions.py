#!/usr/bin/env python
from __future__ import print_function

import os, sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

# Magic spell to make sure Rosetta python libs are on the PYTHONPATH:
#sys.path.append( os.path.abspath( sys.path[0] ) )
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from rosetta_py.io.mdl_molfile import *
#from rosetta_py.utility.rankorder import argmin
#from rosetta_py.utility import r3
from fragment_functions import *

# Features from Python 2.5 that we want to use:
if not hasattr(__builtins__, "any"):
    def any(itr):
        for el in itr:
            if el: return True
        return False

if not hasattr(__builtins__, "all"):
    def all(itr):
        for el in itr:
            if not el: return False
        return True

# Lookup functions for parent residues
def map_parent_resnames_three_to_one( key):
    d = {
        'GLY':'G',
        'CYS':'C',
        'ASP':'D',
        'SER':'S',
        'GLN':'Q',
        'LYS':'K',
        'ILE':'I',
        'THR':'T',
        'PHE':'F',
        'ASN':'N',
        'HIS':'H',
        'LEU':'L',
        'ARG':'R',
        'TRP':'W',
        'ALA':'A',
        'VAL':'V',
        'GLU':'E',
        'TYR':'Y',
        'PRO':'P',
        'MET':'M'
    }
    return d.get( key)


def map_parent_resnames_one_to_three( key):
    d = {

        'G':'GLY',
        'C':'CYS',
        'D':'ASP',
        'S':'SER',
        'Q':'GLN',
        'K':'LYS',
        'I':'ILE',
        'T':'THR',
        'F':'PHE',
        'N':'ASN',
        'H':'HIS',
        'L':'LEU',
        'R':'ARG',
        'W':'TRP',
        'A':'ALA',
        'V':'VAL',
        'E':'GLU',
        'Y':'TYR',
        'P':'PRO',
        'M':'MET'
    }
    return d.get( key)


def write_ligand_kinemage(f, molfile):
    if not isinstance(f, file): f = open(f, 'w')
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
    f.write("@labellist {PDB names} color= white off\n")
    for a in molfile.atoms: f.write("{%s} %.3f %.3f %.3f\n" % (a.pdb_name, a.x, a.y, a.z))
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

def write_param_file(f, molfile, name, frag_id, base_confs, max_confs):
    '''Writes a Molfile object to a file.
    f may be a file name or file handle.
    base_confs is the number of heavy-atom conformations generated by e.g. Omega
    max_confs is the maximum number of conformations desired after adding proton rotation
        The total number of confs may still be larger than max_confs,
        but at least we'll skip -ex# extra sampling of proton chis.
    '''
    close_file = False
    if not isinstance(f, file):
        f = open(f, 'w')
        close_file = True
    if frag_id == 1 and len(name) > 2: name = "%3.3s" % name
    else: name = "%2.2s%1i" % (name, frag_id)
    f.write("NAME %s\n" % name)
    f.write("IO_STRING %3.3s %1i\n" % (name, frag_id))
    f.write("TYPE LIGAND\n")
    f.write("AA UNK\n")
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
        f.write("ATOM %-4s %-4s %-4s %.3f\n" % (atom.name, atom.ros_type, atom.mm_type, atom.partial_charge))
        for a2 in atom.children: write_atoms(a2)
    write_atoms(root_atom)
    for bond in bonds:
        if int(bond.order) == int(4):
            f.write("BOND_TYPE %-4s %-4s %-4s\n" % (bond.a1.pdb_name, bond.a2.pdb_name, 'ARO'))
        else:
            f.write("BOND_TYPE %-4s %-4s %-4s\n" % (bond.a1.pdb_name, bond.a2.pdb_name, bond.order))
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
        print( "WARNING: skipping extra samples for proton chis; would give %i conformers" % num_H_confs )
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
    (nbr, nbr_dist) = choose_neighbor_atom(molfile, frag_id)
    f.write("NBR_ATOM %s\n" % nbr.name)
    f.write("NBR_RADIUS %f\n" % nbr_dist)
    # Convention seems to be a depth-first traversal from the root.
    # I don't know whether this matters, but it's the easy way anyhow.

    # Write formal charge
    formal_charge = int(round(sum(a.partial_charge for a in atoms if not (a.poly_upper or a.poly_lower))))
    if formal_charge != 0:
        if formal_charge < 0:
            f.write("NET_FORMAL_CHARGE %d\n"%formal_charge)
        else:
            f.write("NET_FORMAL_CHARGE +%d\n"%formal_charge)

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
            prev_sibling = conn.a2 # if we allowed multiple connects from one parent
    write_icoords(root_atom)
    #
    # XXX-FIXME:  still might need TYPE, PROPERTIES, FIRST_SIDECHAIN_ATOM, ACT_COORD_ATOMS
    #
    if close_file: f.close()

def get_bonded_distance_map(atoms, bonds, root_atom):
    bonded_distance_map = {}
    bonded_distance_map[ root_atom ] = 0

    def update(ref, new):
        '''Helper function to do the update -- will attempt to find the minimum distance'''
        if ref in bonded_distance_map:
            d = bonded_distance_map[ ref ] + 1
            bonded_distance_map[ new ] = min( d, bonded_distance_map.get(new, 9999999999) )

    for ii in range(len(atoms)): # At worst, we should add one atom per cycle
        for b in bonds:
            update( b.a1, b.a2 )
            update( b.a2, b.a1 )
        if len(bonded_distance_map) == len(atoms):
            break # Early out if we've assigned everything.

    return bonded_distance_map

def write_poly_param_file(f, molfile, name, frag_id, peptoid, parent):
    '''Writes a Molfile object to a file.
    f may be a file name or file handle.
    '''
    close_file = False
    if not hasattr(f, 'write'):
        f = open(f, 'w')
        close_file = True
    else: name = "%2.2s" % (name)
    f.write("NAME %s\n" % name)
    f.write("IO_STRING %3.3s X\n" % (name))
    f.write("TYPE POLYMER\n")
    f.write("AA UNK\n")

    # Specify backbone and sidechain parent rotamers
    if parent is not None:
        if len(parent) == int(3) and map_parent_resnames_three_to_one( parent) is not None:
            f.write("BACKBONE_AA " + str( parent) + "\n")
            f.write("ROTAMER_AA " + str( parent) + "\n")
        elif len(parent) == int(1) and map_parent_resnames_one_to_three( parent) is not None:
            parent_code = map_parent_resnames_one_to_three( parent)
            f.write("BACKBONE_AA " + str( parent_code) + "\n")
            f.write("ROTAMER_AA " + str( parent_code) + "\n")
        else:
            print("WARNING: Invalid parent specified. Check that your canonical amino acid 1- or 3-letter code is correct.")


    atoms = [atom for atom in molfile.atoms if atom.poly_ignore == False]
    bonds = [bond for bond in molfile.bonds if bond.a1.poly_ignore == False and bond.a2.poly_ignore == False]
    root_atoms = [atom for atom in atoms if atom.is_root]
    assert(len(root_atoms) == 1)
    root_atom = root_atoms[0]

    # Reorder atoms based on (bonded) distance to root
    bonded_distance_map = get_bonded_distance_map(atoms, bonds, root_atom)
    atoms.sort( key=lambda a: bonded_distance_map[a] )
    atoms.sort( key=lambda a: a.poly_backbone, reverse=True  ) # Put backbone atoms in front
    atoms.sort( key=lambda a: a.is_H ) # Put hydrogens to the back
    # Reorder bonds based on the atom ordering
    bonds.sort( key=lambda b: max( atoms.index(b.a1), atoms.index(b.a2) ) )

    # write atoms
    for atom in atoms:
        if atom.poly_lower != True and atom.poly_upper != True:
            # Halogenated (Cl, Br) substituents of large sidechains may increase character count; if so, strip leading whitespace
            if len(atom.pdb_name) < int(5):
                f.write("ATOM %-4s %-4s %-4s %.3f\n" % (atom.pdb_name, atom.ros_type, atom.mm_type, atom.partial_charge))
            else:
                f.write("ATOM %-4s %-4s %-4s %.3f\n" % (atom.pdb_name.lstrip(), atom.ros_type, atom.mm_type, atom.partial_charge))
    # write bonds
    for bond in bonds:
        if bond.a1.poly_lower != True and bond.a1.poly_upper != True and bond.a2.poly_lower != True and bond.a2.poly_upper != True:
            if bond.a1.elem == 'H':
                f.write("BOND %-4s %-4s\n" % (bond.a2.pdb_name, bond.a1.pdb_name))
            elif int(bond.order) == int(4):
                # Write aromatic bond types explicitly
                f.write("BOND_TYPE %-4s %-4s %-4s\n" % (bond.a1.pdb_name, bond.a2.pdb_name, 'ARO'))
            else:
                # Is there a reason why bond order is not explicitly written by default?
                f.write("BOND_TYPE %-4s %-4s %-4s\n" % (bond.a1.pdb_name, bond.a2.pdb_name, bond.order))
                #f.write("BOND %-4s %-4s\n" % (bond.a1.pdb_name, bond.a2.pdb_name))

    lower_connects = []
    upper_connects = []
    for bond in bonds:
        if bond.a1.poly_lower: lower_connects.append( bond.a2.pdb_name )
        if bond.a2.poly_lower: lower_connects.append( bond.a1.pdb_name )
        if bond.a1.poly_upper: upper_connects.append( bond.a2.pdb_name )
        if bond.a2.poly_upper: upper_connects.append( bond.a1.pdb_name )

    if len(lower_connects) != 1:
        print("WARNING: One and only one non-ignored atoms should be bonded to the POLY_LOWER atom. Found:", lower_connects)
    if len(upper_connects) != 1:
        print("WARNING: One and only one non-ignored atoms should be bonded to the POLY_UPPER atom. Found:", upper_connects)

    # write upper and lower connect
    f.write("LOWER_CONNECT "+lower_connects[0].strip()+"\n")
    f.write("UPPER_CONNECT "+upper_connects[0].strip()+"\n")

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
            # If we connect backbone to backbone across the bond, then we're not a proper "rotatable" bond
            if (d.poly_backbone == True or c.poly_backbone == True) and (b.poly_backbone == True or a.poly_backbone == True):
                continue
            yield (bond, a, b, c, d)
    def is_sp2_proton(a, b, c, d):
        '''Crude guestimate of H-bonding geometry'''
        # i.e., assume sp2 if the H's grandparent has anything other than single bonds to it.
        return ((a.is_H and len([bnd for bnd in c.heavy_bonds if bnd.order != Bond.SINGLE]) > 0)
             or (d.is_H and len([bnd for bnd in b.heavy_bonds if bnd.order != Bond.SINGLE]) > 0))

    sorted_bonds = list(bonds) # make a copy then sort in place
    sorted_bonds.sort(key=lambda b: b.is_proton_chi, reverse=False)

    # hacky check to see if chis in correct order else revese the order
    if not peptoid:
        rotatable_bonds = list( rot_bond_iter(sorted_bonds) )
        rotatable_bonds.reverse()
        for bond, a, b, c, d in rotatable_bonds:
            if a.poly_upper == False and b.poly_upper == False and c.poly_upper == False and d.poly_upper == False and a.poly_lower == False and b.poly_lower == False and c.poly_lower == False and d.poly_lower == False:
                #print( d.poly_n_bb, d.pdb_name, c.pdb_name, b.pdb_name, a.pdb_name )
                if d.poly_n_bb:
                    sorted_bonds.reverse()
                    print( "REVERSED!!!" )
                    break
                else:
                    break

    num_chis = 0
    for bond, a, b, c, d in rot_bond_iter(sorted_bonds):
        if a.poly_upper == False and b.poly_upper == False and c.poly_upper == False and d.poly_upper == False and \
                a.poly_lower == False and b.poly_lower == False and c.poly_lower == False and d.poly_lower == False:

            # Do not include backbone chi (where "both sides" of the bond are backbone)
            if (d.poly_backbone == True or c.poly_backbone == True) and (b.poly_backbone == True or a.poly_backbone == True):
                continue
            # in peptoid case, there should C backbone should not be in any CHI angles
            if peptoid:
                if d.pdb_name.strip() == "C" or c.pdb_name.strip() == "C" or b.pdb_name.strip() == "C" or a.pdb_name.strip() == "C":
                    continue

            num_chis += 1
            # "Direction" of chi definition appears not to matter,
            # but we follow the convention of amino acids (root to tip)
            f.write("CHI %i %-4s %-4s %-4s %-4s\n" % (num_chis, d.pdb_name, c.pdb_name, b.pdb_name, a.pdb_name))
            if bond.is_proton_chi:
                # Only expand proton chis with extra samples if doing so won't generate tens of thousands of confs.
                extra = "1 20"
                if is_sp2_proton(a, b, c, d):
                    f.write("PROTON_CHI %i SAMPLES 2 0 180 EXTRA %s\n" % (num_chis, extra))
                else:
                    f.write("PROTON_CHI %i SAMPLES 3 60 -60 180 EXTRA %s\n" % (num_chis, extra))

    # determine first side chain atom order of atoms should be n ca c o upper lower [side chain heavys] [hydrogens]
    non_bb_heavy_atoms = [a for a in atoms if a.poly_backbone == False and a.poly_lower == False and a.poly_upper == False and a.poly_ignore == False and a.elem != 'H']
    if len(non_bb_heavy_atoms) > 0:
        first_sidechain_atom = non_bb_heavy_atoms[0]
    else:
        first_sidechain_atom = None

    # choose neighbor atom first CB else CA
    nbr_list = [i for i,a in enumerate(molfile.atoms) if a.pdb_greek_dist == 'B']
    nbr_atom = None
    if len(nbr_list) >= 1:
        nbr_atom_index = nbr_list[0]
        nbr_atom = atoms[nbr_atom_index]
    else:
        for i,a in enumerate(molfile.atoms):
            if a.poly_ca_bb == True:
                nbr_atom_index = i
                nbr_atom = a
                break
    if nbr_atom is None and first_sidechain_atom is not None:
        # We kinda assume that this is attached to the backbone -- perhaps we should be more explicit about it
        nbr_atom = first_sidechain_atom
        nbr_atom_index = [ i for i, a in enumerate(atoms) if a.pdb_name == nbr_atom.name ][ 0 ] # Should just be the one
    if nbr_atom is None:
        print("WARNING: Residue does not have a sidechain -- picking arbitrary neighbor atom from backbone")
        # Attempt to pick an atom in the "middle" of the backbone
        backbone_atoms = [ a for a in atoms if a.poly_backbone and not a.is_H and not a.poly_upper and not a.poly_lower ]
        nbr_atom = backbone_atoms[ len(backbone_atoms) // 2 ]
        nbr_atom_index = [ i for i, a in enumerate(atoms) if a.pdb_name == nbr_atom.name ][ 0 ] # Should just be the one
        # There's probably a better way of dealing with this edge case

    f.write("NBR_ATOM %s\n" % nbr_atom.pdb_name)

    #Using atoms here instead of molfile.atoms to ensure poly_ignore stays ignored
    na = len(atoms) # Number of Atoms
    all_all_dist = [ [1e100] * na for i in range(na) ]
    nbrs = dict([ (a,set()) for a in atoms ])
    for a in atoms:
        nbrs[a].update([b.a2 for b in a.bonds if not b.a2.poly_ignore])
    for i in range(0,na):
        all_all_dist[i] = dijkstra( start = atoms[i], nodes = atoms, nbr = lambda a: nbrs[a], dist = r3.distance )
    nbr_dist = max(all_all_dist[nbr_atom_index])
    f.write("NBR_RADIUS %f\n" % nbr_dist)

    #Write formal charge
    formal_charge = int(round(sum(a.partial_charge for a in atoms if not (a.poly_upper or a.poly_lower))))
    if formal_charge != 0:
        if formal_charge < 0:
            f.write("NET_FORMAL_CHARGE %d\n"%formal_charge)
        else:
            f.write("NET_FORMAL_CHARGE +%d\n"%formal_charge)

    if first_sidechain_atom is not None:
        f.write("FIRST_SIDECHAIN_ATOM %s\n" % first_sidechain_atom.pdb_name)
    else:
        f.write("FIRST_SIDECHAIN_ATOM NONE\n")

    # Specify rama prepro info based on parent
    if parent is not None:
        if len(parent) == int(3) and map_parent_resnames_three_to_one( parent) is not None:
            f.write("RAMA_PREPRO_RESNAME " + str(parent) + "\n")
            f.write("RAMA_PREPRO_FILENAME all.ramaProb prepro.ramaProb" + "\n")
        elif len(parent) == int(1) and map_parent_resnames_one_to_three( parent) is not None:
            parent_code = map_parent_resnames_one_to_three( parent)
            f.write("RAMA_PREPRO_RESNAME " + str(parent_code) + "\n")
            f.write("RAMA_PREPRO_FILENAME all.ramaProb prepro.ramaProb" + "\n")

    # properties
    properties_list = []
    for line in molfile.footer:
        if line.startswith("M  POLY_PROPERTIES"):
            properties_list = line.split()[2:]
    f.write("PROPERTIES")
    for p in properties_list: f.write(" %s" % p)
    f.write("\n")

    # Reporting
    print("WRITING ICOOR_INTERNAL: ")

    def write_icoords(a):
        #print("ICOOR_INTERNAL   %-4s %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s" % (a.pdb_name, a.phi, a.theta, a.d, a.input_stub1.pdb_name, a.input_stub2.pdb_name, a.input_stub3.pdb_name))
        f.write("ICOOR_INTERNAL   %-4s %11.6f %11.6f %11.6f  %-4s  %-4s  %-4s\n" % (a.pdb_name, a.phi, a.theta, a.d, a.input_stub1.pdb_name, a.input_stub2.pdb_name, a.input_stub3.pdb_name));
        for child in a.children:
            if not child.poly_ignore:
                write_icoords(child)

    write_icoords( root_atom )

    # End
    if close_file: f.close()


def write_pdb_rotamers( f, pdb_file_name):
    # Open params file
    close_file = False
    if not hasattr(f, 'write'):
        f = open(f, 'a')
        close_file = True
    f.write( "PDB_ROTAMERS " + str(pdb_file_name))
    # End
    if close_file: f.close()


def write_ligand_pdb(f, molfile_tmpl, molfile_xyz, resname, ctr=None, chain_id='X'):
    '''Writes a PDB file with a series of residues representing one ligand conformation.
    The topology (atom names, fragment divisions, etc) are taken from molfile_tmpl,
    while the actual XYZ coordinates are taken from molfile_xyz.
    resname provides the first two characters of the residue name.
    f may be a file name or file handle.'''
    if not hasattr( f, 'write' ): f = open(f, 'w')
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
        for atom_tmpl in atoms:
            # skip this atom if it is ignored in polymer residue case
            if atom_tmpl.poly_ignore or atom_tmpl.poly_lower or atom_tmpl.poly_upper:
                continue
            atom_xyz = molfile_xyz.atoms[ ai[atom_tmpl]-1 ]
            xyz = r3.add(atom_xyz, ctr)
            atom_num += 1
            if len(frag_ids) == 1 and len(resname) > 2 and atom_tmpl.poly_ignore == False:
                f.write("HETATM%5i %-4.4s %3.3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  \n"
                    % (atom_num, atom_tmpl.pdb_name, resname,          chain_id, frag_id, xyz.x, xyz.y, xyz.z, 1.0, 20.0, atom_tmpl.elem))
            else:
                f.write("HETATM%5i %-4.4s %2.2s%1i %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  \n"
                    % (atom_num, atom_tmpl.pdb_name, resname, frag_id, chain_id, frag_id, xyz.x, xyz.y, xyz.z, 1.0, 20.0, atom_tmpl.elem))
    f.write("TER"+(" "*77)+"\n")
    f.close()

