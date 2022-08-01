#!/usr/bin/env python
from __future__ import print_function

import os, sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

from optparse import OptionParser

# Magic spell to make sure Rosetta python libs are on the PYTHONPATH:
#sys.path.append( os.path.abspath( sys.path[0] ) )
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from rosetta_py.io.mdl_molfile import *
#from rosetta_py.utility.rankorder import argmin
#from rosetta_py.utility import r3

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

def add_fields_to_bonds(bonds):
    '''Adds a bunch of member variable placeholders that we use.'''
    for bond in bonds:
        bond.can_rotate = False     # true for single bonds not in rings
        bond.is_proton_chi = False  # true for bonds that rotate OH hydrogens, etc
        bond.connection_id = 0      # non-zero id if bond crosses fragments
        # Remember we have to update mirror too!
        bond.mirror.can_rotate      = bond.can_rotate
        bond.mirror.is_proton_chi   = bond.is_proton_chi
        bond.mirror.connection_id   = bond.connection_id
        bond.poly_ignore = False    # convience boolean

def check_bond_count(atoms):
    '''Safety check for structures with stupid numbers of bonds to atoms, at Florian's request.'''
    def valence(bond):
        if bond.order == Bond.AROMATIC: return 1.5
        else: return float(bond.order)
    for atom in atoms:
        if atom.is_H and len(atom.bonds) > 1:
            raise ValueError("Atom %s is a hydrogen with >1 bonds" % atom.name)
        # Valence 4.5 for e.g. carbon comes up in peptide bonds and at the joins in multi-ring systems.
        if sum([valence(bond) for bond in atom.bonds]) > 4.5:
            print( "WARNING:  atom %s has valence > 4" % atom.name )

def check_aromaticity(bonds):
    '''Safety check for Kekule structures (alternating single/double bonds)
    rather than bonds described explicitly as aromatic.'''
    aro_bonds = len([b for b in bonds if b.order == Bond.AROMATIC])
    dbl_bonds = len([b for b in bonds if b.order == Bond.DOUBLE])
    if aro_bonds == 0 and dbl_bonds > 0:
        print( "WARNING:  structure contains double bonds but no aromatic bonds" )
        print( "  Aromatic bonds must be identified explicitly --" )
        print( "  alternating single/double bonds (Kekule structure) won't cut it." )
        print( "  This warning does not apply to you if your molecule really isn't aromatic." )

def assign_rotatable_bonds(bonds):
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

