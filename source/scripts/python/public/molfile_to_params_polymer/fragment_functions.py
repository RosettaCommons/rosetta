#!/usr/bin/env python
from __future__ import print_function

import os, sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

# Magic spell to make sure Rosetta python libs are on the PYTHONPATH:
#sys.path.append( os.path.abspath( sys.path[0] ) )
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

#from rosetta_py.io.mdl_molfile import *
from rosetta_py.utility.rankorder import argmin
from rosetta_py.utility import r3

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

def atom_num( desig, m ):
    try:
        return int(desig) - 1
    except ValueError:
        # Atom name?
        for ii, a in enumerate(m.atoms):
            if desig == a.name:
                return ii

    raise ValueError("Cannot find atom designation {}, should be integer or one of {}.".format( desig, [a.name for a in m.atoms] ) )

def assign_rigid_ids(atoms):
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

def fragment_ligand(molfile):
    '''Sets Atom.fragment_id, Atom.conn_bonds, and Bond.connection_id.
    Returns number of fragments created, i.e. the largest valid fragment id.'''
    remaining_bonds = set(molfile.bonds) # copy
    # Delete all split bonds from remaining_bonds
    # and number the connections they leave behind
    num_conn_id = 0
    for line in molfile.footer:
        if not line.startswith("M SPLT"): continue
        fields = line.split()
        atom1 = molfile.atoms[ atom_num(fields[2], molfile) ]
        atom2 = molfile.atoms[ atom_num(fields[3], molfile) ]
        bond_to_remove = None
        for b in remaining_bonds:
            if((b.a1 == atom1 and b.a2 == atom2)
            or (b.a1 == atom2 and b.a2 == atom1)):
                bond_to_remove = b
                break
        if b is None:
            raise ValueError("Cannot find bond to split between %s and %s" % (atom1.name, atom2.name))
        #elif b.can_rotate:
            #raise ValueError("Shouldn't split ROTATABLE bond between %s and %s" % (atom1.name, atom2.name))
        else:
            if b.can_rotate:
                print( "WARNING: spliting ROTATABLE bond between %s and %s" % (atom1.name, atom2.name) )
            print( "Split bond between %s and %s" % (atom1.name, atom2.name) )
            num_conn_id += 1
            bond_to_remove.connection_id = num_conn_id
            bond_to_remove.mirror.connection_id = bond_to_remove.connection_id
            bond_to_remove.a1.conn_bonds.append(bond_to_remove)
            bond_to_remove.a2.conn_bonds.append(bond_to_remove.mirror)
            remaining_bonds.remove(bond_to_remove)
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
    #    print( atom.name, atom.fragment_id )
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
            raise ValueError("Cannot create more than one connection at same atom (%s)" % atom.name)
        # Only check one permutation b/c the other gets caught going backwards along the same bond.
        for conn in atom.conn_bonds:
            pair1 = (conn.a1.fragment_id, conn.a2.fragment_id)
            #pair2 = (conn.a2.fragment_id, conn.a1.fragment_id)
            if pair1 in frag_frag_conns:
                #raise ValueError("Cannot create multiple connections between fragments %i and %i" % pair1)
                print( "WARNING: Multiple connections between fragments (%i and %i) **NOT CURRENTLY SUPPORTED** by Rosetta!" % pair1 )
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
            print( "Fragment %i: %s" % (frag_id, [a.name for a in frag_atoms]) )
            raise ValueError("Fragment %i has %i atoms; merge with another fragment or add virtual atoms to make 3 total" % (frag_id, num_atoms))
        if not (7 <= num_atoms <= 24):
            print( "WARNING: fragment %i has %i total atoms including H; protein residues have 7 - 24 (DNA: 33)" % (frag_id, num_atoms) )
        if not (4 <= num_heavy_atoms <= 22):
            print( "WARNING: fragment %i has %i non-H atoms; protein residues have 4 - 14 (DNA: 22)" % (frag_id, num_heavy_atoms) )
        if num_rot_bonds > 4:
            print( "WARNING: fragment %i has %i rotatable bonds; protein residues have 0 - 4" % (frag_id, num_rot_bonds) )
    print( "Average %.1f atoms (%.1f non-H atoms) per fragment" % (
        float(len(molfile.atoms)) / float(num_frag_id), float(len([a for a in molfile.atoms if not a.is_H])) / float(num_frag_id)) )
    # Average stats tabulated by IWD from Richardson's Top500 database
    print( "(Proteins average 15.5 atoms (7.8 non-H atoms) per residue)" )
    return num_frag_id

def build_fragment_trees(molfile):
    '''Assigns a root atom for each fragment and parents and children.'''
    # Assign root atoms based on instructions in the molfile
    for line in molfile.footer:
        # Standard MDL style is with a space, but KWK has used "MROOT" in the past. Babel uses 2 spaces.
        if   line.startswith("M  ROOT"): molfile.atoms[ atom_num(line.split()[2], molfile) ].is_root = True
        elif line.startswith("M ROOT"):  molfile.atoms[ atom_num(line.split()[2], molfile) ].is_root = True
        elif line.startswith("MROOT"):   molfile.atoms[ atom_num(line.split()[1], molfile) ].is_root = True
    for frag_id in set([a.fragment_id for a in molfile.atoms]):
        # If we want to have a default way of choosing the root atom, this is the place:
        root_atoms = [a for a in molfile.atoms if a.fragment_id == frag_id and a.is_root]
        if len(root_atoms) == 0:
            print( "WARNING:  no root atom specified, using auto-selected NBR atom instead." )
            (nbr, nbr_dist) = choose_neighbor_atom(molfile, frag_id)
            nbr.is_root = True
            root_atoms = [nbr]
        elif len(root_atoms) != 1:
            raise ValueError("You must have no more than one 'M ROOT' record in the molfile per ligand fragment")
        root_atom = root_atoms[0]
        # Protein residues appear to go depth first, so that all chi angles ride on each other.
        # Depth first assignment -- leads to very deep trees
        def tree_dfs(parent):
            def get_atom_num(elem):
                elem_atom_num = {'H':1,'B':5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'NA': 11, 'MG': 12, 'P':15, 'S':16, 'CL':17, 'K':19, 'CA':20,
                                 'FE':26, 'ZN':30, 'BR':35, 'I':53, 'SE':34, 'X':-1}
                try:
                    return elem_atom_num[elem]
                except:
                    print("WARNING: element type %s not recognized, atom tree construction may be incorrect."%elem)
                    return 99
            # Want to visit non-H children first
            tmp_children = [b.a2 for b in parent.bonds]
            def sort_key(a):
                # Priority to higher valued items (reverse sort)
                return (not a.poly_ignore, # Priority to non-ignore
                        not ( a.poly_upper or a.poly_lower ), # Priority to non-connect
                        a.poly_backbone, # Priority to backbone
                        get_atom_num(a.elem), # Priority to higher atom number
                        max([get_atom_num(b.a2.elem) for b in a.bonds]), # Priority to bonded to heavier MW
                        len([b for b in a.bonds if not b.a2.is_H]), ) # Number of non-hydrogen bonds
            tmp_children.sort(key=sort_key, reverse=True)

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

def handle_polymer_special_cases(child, molfile):
    # There should be a better way of doing this
    nbb = cabb = cbb = obb = upper = lower = None
    for a in molfile.atoms:
        if a.poly_n_bb == True: nbb = a
        if a.poly_ca_bb == True: cabb = a
        if a.poly_c_bb == True: cbb = a
        if a.poly_o_bb == True: obb = a
        if a.poly_upper == True: upper = a
        if a.poly_lower == True: lower = a

    # For alpha amino acids
    if child.poly_n_bb:
        child.input_stub1, child.input_stub2, child.input_stub3 = nbb, cabb, cbb
    elif child.poly_ca_bb:
        child.input_stub1, child.input_stub2, child.input_stub3 = nbb, cabb, cbb
    elif child.poly_c_bb:
        child.input_stub1, child.input_stub2, child.input_stub3 = cabb, nbb, cbb
    elif child.poly_o_bb:
        child.input_stub1, child.input_stub2, child.input_stub3 = cbb, cabb, upper
    elif child.poly_upper and None not in (cbb, cabb, nbb):
        child.input_stub1, child.input_stub2, child.input_stub3 = cbb, cabb, nbb
    elif child.poly_lower and None not in (nbb, cabb, cbb):
        child.input_stub1, child.input_stub2, child.input_stub3 = nbb, cabb, cbb
    else:
        return

def assign_internal_coords(molfile):
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
            #print( "stubs", [x.name for x in (me, me.stub1, me.stub2, me.stub3)] )
            # Assign input stubs to children and calculate internal coords
            prev_sibling = me.stub3
            parent = me # rename to make logic clearer
            # Have to store input_stub atoms b/c they're written to params file

            # Descend into backbone atoms first, then sidechain
            for child in sorted( parent.children, key=lambda a: a.poly_backbone, reverse=True ):
                child.input_stub1 = parent.stub1
                child.input_stub2 = parent.stub2
                child.input_stub3 = prev_sibling # for first child, this is parent.stub3
                # Special case for second child of the root
                if parent.is_root and prev_sibling == parent.stub2:
                    #print( "activate second child case! stub3 =", parent.stub3.name )
                    child.input_stub3 = parent.stub3

                handle_polymer_special_cases(child, molfile)

                #print( "input_stubs", [x.name for x in (child, child.parent, child.input_stub1, child.input_stub2, child.input_stub3)] )
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
        handle_polymer_special_cases(root_atom, molfile) # In case root is a special case (which it probably is)
        # Root has dummy values for d, theta, phi
        root_atom.d     = 0.0
        root_atom.theta = 0.0
        root_atom.phi   = 0.0
    # end loop over all fragments

def calc_internal_coords(child, input_stub1, input_stub2, input_stub3):
    '''Returns (d, theta, phi) for a point given it's three input stub points.'''
    #print( "calc_internal_coords", [x.name for x in (child, input_stub1, input_stub2, input_stub3)] )
    # Now actually calculate spherical internal coordinates
    # There's some weird "flip_stub" logic for first child of the root
    # iff theta is to be kept fixed (BondedAtom.cc) but I don't think it matters here.
    #
    # parent == parent.stub1 == child.input_stub1, always
    # (except for CONNECT atoms, where we use different stubs!!)
    d = r3.distance(child, input_stub1)
    if d < 1e-2: # very small d;  theta, phi don't make sense
        print( "WARNING: very small d=%f for %s" % (d, child.name) )
        theta = 0.0
        phi = 0.0
    else:
        theta = r3.angle( r3.from_to(input_stub2,input_stub1), r3.from_to(input_stub1,child) )
        if theta < 1e-2 or theta > 180 - 1e-2:
            # This always happens for first child of root:
            #print( "WARNING: nearly parallel theta=%f for %s" % (theta, child.name) )
            phi = 0.0
        else:
            #print("PHI: " + child.pdb_name + " " + input_stub1.pdb_name + " " + input_stub2.pdb_name + " " + input_stub3.pdb_name)
            phi = r3.dihedral(child, input_stub1, input_stub2, input_stub3)
    return (d, theta, phi)

def choose_neighbor_atom(molfile, frag_id):
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
        for i in range(0,na):
            all_all_dist[i] = dijkstra(
                start = atoms[i],
                nodes = atoms,
                nbr = lambda a: nbrs[a],
                dist = r3.distance
            )
    # The "best" neighbor atom is the one with
    # the smallest {max distance to all other atoms}.
    max_dist = [ max(row) for row in all_all_dist ]
    best_idx = argmin(max_dist)
    return (atoms[best_idx], max_dist[best_idx])

def dijkstra(start, nodes, nbr, dist):
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
    #print( "queue is ",queue #DEBUG )
    # Allows lookup of best distance by name
    shortest = dict([ (q[NODE],q) for q in queue ])
    #print( "start shortest is ", shortest[start][DIST] ) #DEBUG
    shortest[start][DIST] = 0 # start
    while len(queue) > 0:
        # return the index of the smallest element of the queue
        curr_idx = argmin(queue) # on first pass this is start
        curr = queue.pop(curr_idx)[NODE] # pop out the node with smallest
        curr_shortest = shortest[curr][DIST]
        for n in nbr(curr):
            new_dist = curr_shortest + dist(curr,n)
            if new_dist < shortest[n][DIST]:
                shortest[n][DIST] = new_dist
    return [ shortest[node][DIST] for node in nodes ]

