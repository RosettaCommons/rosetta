#!/usr/bin/env python3
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
'''
Contains helper functions that defines Rosetta-style topology for given Molecule including:
- AtomTree, Chi angles, Ring, Bond order, and so on.

Author: Hahnbeom Park and Frank DiMaio 
'''

import sys
import math
import numpy
import numpy.linalg
import scipy.sparse.csgraph
import scipy.linalg
from Types import *
from BasicClasses import AtomClass, RingClass
from AtomTypeClassifier import AtomTypeClassifier
from utils import dihedral

def setup(mol,option): # mol: Molecule type
    '''Main function setting up all the toplogy-related things of the input Molecule(=="mol") '''
    mol.vatms = []
    mol.angles = []
    mol.torsions = []
    mol.nbratom = -1
    mol.nbrradius = 0.0
    mol.atms_aro = []
    mol.atms_puckering = []
    mol.atms_ring = [] #including aromatic ring
    mol.ATorder = []
    mol.rings_aro = [] #only aromatic
    mol.rings_pucker = [] #puckering
    mol.rings_sugar = [] #subset of puckering rings
    mol.rings = [] #aromatic + puckering

    assign_bonds(mol)
    assign_angles(mol)
    assign_torsions(mol)
    detect_rings(mol,option)

    assign_hybridization_if_missing(mol) # should come after ring assignment
    setup_nbratom(mol) #should come after ring assignments
    define_icoord(mol) #should come after nbr setup

    classifier = AtomTypeClassifier()
    classifier.apply_to_molecule(mol)
    classifier.assert_H(mol)

    # should come after atom types classified
    define_conjugation(mol)
    define_rotable_torsions(mol)

# Functions from here are pieces of setup

def setup_nbratom(mol):
    '''Setup residue nbratom & nbrradius field'''
    com = numpy.sum(mol.xyz,axis=0) / len(mol.atms)
    dists = numpy.sum(numpy.square(mol.xyz - com[None,:]),axis=-1) # dist_squared
    maxdis = math.sqrt(max(dists)) 

    # 1) dont let puckering ring be root
    dists[numpy.in1d(numpy.arange(len(mol.atms)), mol.atms_puckering)] = 99999.0

    nbonds = [len(atm.bonds) for atm in mol.atms]
    for ring in mol.rings: #correct for cut_bonds
        (i,j) = ring.cut_bond
        if i == None: continue
        nbonds[i] -=1
        nbonds[j] -=1

    # 3) also exclude hydrogens
    # prioritize ring-atom with >2 bonds
    for i,atm in enumerate(mol.atms):
        if nbonds[i] <=1: dists[i] = 99999.0
        elif nbonds[i] ==2: dists[i] = 999.0
        if atm.is_H: dists[i] = 99999.0

    mol.nbratom = numpy.argmin(dists)
    mol.nbrradius = (maxdis+1.5)*2 #safe

def assign_bonds(mol):
    for bond in mol.bonds:
        is_H1 = (mol.atms[bond.atm1].atype==2) 
        is_H2 = (mol.atms[bond.atm2].atype==2) 
        mol.atms[bond.atm1].add_bond(bond.atm2,bond.order)
        mol.atms[bond.atm2].add_bond(bond.atm1,bond.order)

        if is_H1:
            mol.atms[bond.atm2].has_H = True
        if is_H2:
            mol.atms[bond.atm1].has_H = True

        if mol.atms[bond.atm1].atype == 1 and mol.atms[bond.atm2].atype in POLAR_ATOMS:
            mol.atms[bond.atm1].connected_to_polar = True
        if mol.atms[bond.atm2].atype == 1 and mol.atms[bond.atm1].atype in POLAR_ATOMS:
            mol.atms[bond.atm2].connected_to_polar = True

    # sanity check if any atm is not connected to anything
    atms_unconnected = []
    for atm in mol.atms:
        if len(atm.bonds) == 0:
            atms_unconnected.append(atm.name)
    if atms_unconnected != []:
        sys.exit('ERROR: Atom found not connected to any other atom:'+' '.join(atms_unconnected))
        
# Currently used only for 'nh' case: See "0" in SPECIAL_HYBRIDS to see which atypes in mol2 comes here
def assign_hybridization_if_missing(mol):
    n_to_assign = 0
    to_assign = []
    for iatm,atm in enumerate(mol.atms):
        if ATYPES[atm.atype] in ATYPES_HYBRID:
            n_to_assign += 1
            if atm.hyb == 0:
                to_assign.append(iatm)
            
    for iatm in to_assign:
        atm = mol.atms[iatm]
        attached_to_aro = False
        attached_to_nonaro = False
        nH_attached = 0
        for jatm,order in atm.bonds:
            if jatm in mol.atms_aro:
                attached_to_aro = True
            elif mol.atms[jatm].is_H:
                nH_attached += 1
            else:
                attached_to_nonaro = True
        
        if atm.atype == ATYPES.index('N') and attached_to_aro:
            if nH_attached == 2: #aro-NH2
                atm.hyb = 2 
            elif not attached_to_nonaro:
                atm.hyb = 2
            else:
                atm.hyb = 3
        else:
            if len(atm.bonds) < 3:
                atm.hyb = 2
            else:
                #measure improper torsion
                neigh = [j for j,order in atm.bonds][:3]
                imp = dihedral(mol.xyz[neigh[0]],mol.xyz[neigh[1]],mol.xyz[neigh[2]],mol.xyz[iatm])
                if abs(imp) < 10.0: #dev less than 20 degree
                    atm.hyb = 2
                else:
                    atm.hyb = 3
           
        # BELOW FAILS for GAFFtype nh!
        # Use user input structure instead above (there is no other easy way)
        #bondorders = [bond[1] for bond in atm.bonds]
        #maxorder = max(bondorders)
        #if 1 in bondorders:
        #    atm.hyb = 3
        #elif maxorder == 9:
        #    atm.hyb = 9
        #elif maxorder == 2:
        #    atm.hyb = 2
        #elif maxorder == 3:
        #    atm.hyb = 1
        #else:
        #    print('WARNING: Unknown atype in input mol2 for %s; assigning hyb=2'%(atm.name))
        print('Hybridization state for %d:%s assigned as %s'%(iatm,atm.name,atm.hyb))
            
def assign_angles(mol):
    '''Assigns all bond-angles (not only in AtomTree)'''
    for i,bond1 in enumerate(mol.bonds[:-1]):
        for bond2 in mol.bonds[i+1:]:
            if bond1.atm1 == bond2.atm1:
                angle = (bond1.atm2,bond1.atm1,bond2.atm2)
            elif bond1.atm1 == bond2.atm2:
                angle = (bond1.atm2,bond1.atm1,bond2.atm1)
            elif bond1.atm2 == bond2.atm1:
                angle = (bond1.atm1,bond1.atm2,bond2.atm2)
            elif bond1.atm2 == bond2.atm2:
                angle = (bond1.atm1,bond1.atm2,bond2.atm1)
            else:
                continue
            mol.angles.append(angle)

def assign_torsions(mol):
    '''Assigns all torsion angles (not only in AtomTree)'''
    bonds = [(bond.atm1,bond.atm2) for bond in mol.bonds]

    torsions_heavy = []
    torsions_H = []
    orders_heavy = []
    orders_H = []
    if mol.option.verbose:
        print('bonds: ', [(bond.atm1,bond.atm2) for bond in mol.bonds])
    for i,bond1 in enumerate(mol.bonds[:-1]):
        for bond2 in mol.bonds[i+1:]:
            if (bond1.atm1 == bond2.atm1 or bond1.atm1 == bond2.atm2 or \
                bond1.atm2 == bond2.atm1 or bond1.atm2 == bond2.atm2): continue

            if (bond1.atm1,bond2.atm1) in bonds:
                a,b,c,d = (bond1.atm2,bond1.atm1,bond2.atm1,bond2.atm2)
            elif (bond1.atm1,bond2.atm2) in bonds:
                a,b,c,d = (bond1.atm2,bond1.atm1,bond2.atm2,bond2.atm1)
            elif (bond1.atm2,bond2.atm1) in bonds:
                a,b,c,d = (bond1.atm1,bond1.atm2,bond2.atm1,bond2.atm2)
            elif (bond1.atm2,bond2.atm2) in bonds:
                a,b,c,d = (bond1.atm1,bond1.atm2,bond2.atm2,bond2.atm1)
            elif (bond2.atm1,bond1.atm1) in bonds:
                a,b,c,d = (bond2.atm2,bond2.atm1,bond1.atm1,bond1.atm2)
            elif (bond2.atm2,bond1.atm1) in bonds:
                a,b,c,d = (bond2.atm1,bond2.atm2,bond1.atm1,bond1.atm2)
            elif (bond2.atm1,bond1.atm2) in bonds:
                a,b,c,d = (bond2.atm2,bond2.atm1,bond1.atm2,bond1.atm1)
            elif (bond2.atm2,bond1.atm2) in bonds:
                a,b,c,d = (bond2.atm1,bond2.atm2,bond1.atm2,bond1.atm1)
            else:
                continue

            ibond = bonds.index((b,c))
            ## bug fix Sep25 2018: if having only hydrogens at either end
            if mol.atms[a].is_H or mol.atms[d].is_H:
                torsions_H.append((a,b,c,d))
                orders_H.append(mol.bonds[ibond].order)
            else:
                torsions_heavy.append((a,b,c,d))
                orders_heavy.append(mol.bonds[ibond].order)

    # sort in heavy-atom torsions followed by H-only torsions
    #mol.torsion_orderso = orders_heavy + orders_H #unused
    mol.torsions = torsions_heavy + torsions_H

    if mol.option.verbose:
        print('Torsions: ', mol.torsions)

def search_special_biaryl_ring(mol,ring):
    '''Add extra rotatable bonds between ring & nonring functional groups
    (by going through hard-coded cases)'''
    
    biaryl_pivot_extra = []
    # [(A,B),..] : ring=A becomes special biaryl if any of ring=A-B connection is detected
    # Nin definition changed as of Oct242018; All the non-aromatic 1H sp2N goes to NG21
    # non-biaryl (== conjugated): NG21=O*
    special_biaryl_to_ring = [(ACLASS_ID.index('CDp'),ACLASS_ID.index('OG2')),
                              (ACLASS_ID.index('CDp'),ACLASS_ID.index('Oal')),
                              (ACLASS_ID.index('CDp'),ACLASS_ID.index('Oad')),
                              
                              # ring-N=C; (RingN-C)-NG21 are excluded outside here
                              (ACLASS_ID.index('NG21'),ACLASS_ID.index('CR')),
                              (ACLASS_ID.index('NG21'),ACLASS_ID.index('CRp')),
                              (ACLASS_ID.index('NG21'),ACLASS_ID.index('CD')),
                              (ACLASS_ID.index('NG21'),ACLASS_ID.index('CDp')),
                              (ACLASS_ID.index('NG21'),ACLASS_ID.index('CSp')),
                              (ACLASS_ID.index('NG21'),ACLASS_ID.index('CS')),
                              (ACLASS_ID.index('NG21'),ACLASS_ID.index('CS1')),
                              (ACLASS_ID.index('NG21'),ACLASS_ID.index('CS2')),

                              # ring-N-SOn
                              #(ACLASS_ID.index('Nin'),ACLASS_ID.index('SG5')),
                              (ACLASS_ID.index('NG21'),ACLASS_ID.index('SG5')),
                              
                              # ring-amide
                              (ACLASS_ID.index('Nad'),ACLASS_ID.index('CDp')), # treat energy separately

                              # ring-C=N*
                              (ACLASS_ID.index('CDp'),ACLASS_ID.index('Nad')),
                              (ACLASS_ID.index('CDp'),ACLASS_ID.index('Nam2')),
                              (ACLASS_ID.index('CDp'),ACLASS_ID.index('NG2')),
                              (ACLASS_ID.index('CDp'),ACLASS_ID.index('Nin')),

                              # R-CD1-CD1-
                              (ACLASS_ID.index('CD1'),ACLASS_ID.index('CD1')),
                              (ACLASS_ID.index('CD1'),ACLASS_ID.index('CD')),
                              (ACLASS_ID.index('CD'),ACLASS_ID.index('CD')),
                              (ACLASS_ID.index('CD'),ACLASS_ID.index('CD1')),
                              
                              # R-CDx-R
                              (ACLASS_ID.index('CD'),ACLASS_ID.index('CR')),
                              (ACLASS_ID.index('CD1'),ACLASS_ID.index('CR')),
                              (ACLASS_ID.index('CDp'),ACLASS_ID.index('CR')),

                              # ring-C=halogen 
                              (ACLASS_ID.index('CD'),ACLASS_ID.index('F')),
                              (ACLASS_ID.index('CD'),ACLASS_ID.index('Cl')),
                              (ACLASS_ID.index('CD'),ACLASS_ID.index('Br')),
                              (ACLASS_ID.index('CD'),ACLASS_ID.index('I')),
                              ]
    

    for a1 in ring.atms:
        a2_is_special_biaryl = False
        for a2,order in mol.atms[a1].bonds:
            a2class = mol.atms[a2].aclass
            # avoid from ring itself
            #if a2 in mol.atms_aro or a2 in mol.atms_puckering: continue
            if mol.in_same_ring(a1,a2): continue

            # Oct17: skip if not defined as single bond as input!
            if order > 1: continue
            
            further_connected = False
            for a3,order in mol.atms[a2].bonds:
                if a3 == a1: continue
                a3class = mol.atms[a3].aclass
                
                if (a2class,a3class) in special_biaryl_to_ring:
                    a2_is_special_biaryl = True
                    #break
                if len(mol.atms[a3].bonds) > 1:
                    further_connected = True

            if a2_is_special_biaryl and further_connected:
                biaryl_pivot_extra.append((a1,a2))
            
    return biaryl_pivot_extra

def check_duplication(rings,ring_i):
    for ring_j in rings:
        n = 0
        for j in ring_j.atms:
            if j in ring_i: n+=1
        if n == ring_j.natms:
            return ring_j.atms
    return False

def get_path(bond_tree,edge,alt_cut=[]):
    '''helper function in ring detection'''
    i,j = edge
    curr = j
    ring_i = [curr]
    bond_tree[i,j] = bond_tree[j,i] = 0
    if alt_cut != []: bond_tree[alt_cut[0],alt_cut[1]] = bond_tree[alt_cut[1],alt_cut[0]] = 0
    ds, preds = scipy.sparse.csgraph.shortest_path(bond_tree,return_predecessors=True)
    bond_tree[i,j] = bond_tree[j,i] = 1
    if alt_cut != []: bond_tree[alt_cut[0],alt_cut[1]] = bond_tree[alt_cut[1],alt_cut[0]] = 1

    if ds[i,j] > 9999: return False
    
    while (curr != i):
        curr = preds[i,curr]
        ring_i.append(curr)
    return ring_i

def detect_rings(mol,option):
    '''Detects ring & ring types (i.e. aro vs nonaro)'''
    mol.rings = []
    mol.atms_aro = []
    mol.atms_puckering = []
    
    bond_tree = numpy.zeros((len(mol.atms),len(mol.atms)))
    for bond in mol.bonds:
        bond_tree[bond.atm1,bond.atm2] = bond_tree[bond.atm2,bond.atm1] = 1

    mol.tree = scipy.sparse.csgraph.minimum_spanning_tree(bond_tree)
    treesymm = numpy.maximum(mol.tree.toarray(), mol.tree.T.toarray())
    cycle_edges = numpy.transpose(numpy.nonzero(numpy.triu(bond_tree != treesymm)))

    for i,j in cycle_edges:
        ring_i = get_path(bond_tree,(i,j))

        # this sometimes happens when cycle_edge connects two rings
        dupl_ring = check_duplication(mol.rings,ring_i)
        
        if dupl_ring: #search alternative path by cutting an edge of duplicating ring 
            ring_i = get_path(bond_tree,(i,j),alt_cut=(dupl_ring[0],dupl_ring[-1]))
            if ring_i:
                dupl_ring = check_duplication(mol.rings,ring_i)
                if not dupl_ring:
                    mol.rings.append(RingClass(ring_i,cut=(ring_i[-1],ring_i[0])))
                    
        else: #One cycle edge can form alternative ring
            mol.rings.append(RingClass(ring_i,cut=(ring_i[-1],ring_i[0])))
            ring_ii = get_path(bond_tree,(i,j),alt_cut = (ring_i[0],ring_i[1]))
            if ring_ii:
                dupl_ring = check_duplication(mol.rings,ring_ii)
                if not dupl_ring:
                    mol.rings.append(RingClass(ring_ii,cut=(ring_ii[-1],ring_ii[0])))

    for ring in mol.rings:
        classify_ring_type(mol,ring,option)

    mol.atms_ring = mol.atms_aro + mol.atms_puckering

def classify_ring_type(mol,ring,option):
    '''Classifies ring type'''
    is_aro = True

    #for long ring take into account of input geometry
    if ring.natms > 6: 
        ring_xs = [mol.xyz[i] for i in ring.atms]
        ring_xs -= numpy.mean(ring_xs, axis = 0) 
        cov = numpy.cov(ring_xs, rowvar = False)
        evals,_ = scipy.linalg.eigh(cov)
        nonplanarity = numpy.min(evals)
        if (nonplanarity>1e-2):  # be a bit permissive
            is_aro = False
            #print("nonplanarity? ", nonplanarity)
    else: #3~6 membered rings
        nsp2N = 0
        nsp3OS = 0
        nsp3 = 0
        for atm in ring.atms:
            if mol.atms[atm].hyb == 3:
                nsp3 += 1
                for j,order in mol.atms[atm].bonds:
                    if order not in [2,4]: 
                        is_aro = False
                if ATYPES[mol.atms[atm].atype] in ['O','S']:
                    nsp3OS += 1
            elif mol.atms[atm].hyb == 2 and ATYPES[mol.atms[atm].atype] == 'N':
                nsp2N += 1

            # allow exception: 2-nitrogen 5-membered ring or one sp3 O/S ring
            if len(ring.atms) == 5:
                if (nsp3 == nsp3OS) or (nsp2N == 2): 
                    is_aro = True
                else:
                    is_aro = False

    if is_aro:
        ring.type = 4
        mol.rings_aro.append(ring)
        for atm in ring.atms:
            if atm not in mol.atms_aro:
                mol.atms_aro.append(atm)
                        
    else: #otherwise puckering
        ring.type = 2 #non-ringsampling 
        natms = len(ring.atms)
        has_nonsp3 = False
        for atm in ring.atms:
            if mol.atms[atm].hyb != 3:
                has_nonsp3 = True
                break

        if option.opt.report_puckering_chi:
            if (natms > option.opt.longest_puckering_ring):
                ring.type = 3 #long
            elif (natms > 4) and (not option.opt.ring_sampling_sp3_only or not has_nonsp3 ):
                ring.type = 1 #ring-sampling
        
        #sugar
        if len(ring.atms) == 5:
            n_sp3C = 0
            n_O = 0
            for iatm in ring.atms:
                atm = mol.atms[iatm]
                if ATYPES[atm.atype] == 'C' and atm.hyb == 3:
                    n_sp3C += 1
                if ATYPES[atm.atype] == 'O': n_O += 1
            if n_sp3C == 4 and n_O == 1:
                #ring.type = 5
                print('Sugar ring %s: '%mol.mol2file+'  %3s'*5%tuple([mol.atms[iatm].name for iatm in ring.atms]))
        mol.rings_pucker.append(ring)
        
        # puckering also addes into rings (for the biaryl assignment...)
        for atm in ring.atms:
            if atm not in mol.atms_puckering:
                mol.atms_puckering.append(atm)

                
# Scipy-version
def define_icoord(mol):
    '''AtomTree setup using scipy graph construct'''
    nodes,parents = scipy.sparse.csgraph.breadth_first_order(mol.tree, mol.nbratom, directed=False)
    first_children = numpy.zeros_like(parents)
    mol.ATorder = nodes
    
    for i in nodes:
        par_i = parents[i]
        gp_i = -9999 if (par_i==-9999) else parents[par_i]
        ggp_i = -9999 if (gp_i==-9999) else parents[gp_i]

        # 1: root corrections
        if (i == nodes[0] or i == nodes[1] or i == nodes[2]):
            par_i = nodes[0]
            gp_i = nodes[1]
            ggp_i = nodes[2]

        # 2: near-root corrections
        if (gp_i == -9999):
            # 3+ child of root
            gp_i = nodes[1]
            ggp_i = nodes[2]

        if (ggp_i == -9999):
            if (par_i == nodes[1]):
                ggp_i = nodes[2]
            else:
                ggp_i = nodes[1]

        mol.atms[i].root = par_i
        mol.atms[i].groot = (gp_i,ggp_i)

        # ring virtuals
        if (i in mol.atms_puckering):
            for rnum,ring in enumerate(mol.rings):
                (k,l) = ring.cut_bond
                if not (k==i or l==i): continue
                if ring.type != 1: continue
                
                if (i==k):
                    vrt_i = l
                    vtag = "V%dl"%rnum
                else:
                    vrt_i = k
                    vtag = "V%du"%rnum
                    
                # define virtual atoms
                atm = AtomClass(vtag,"X",0,0.0)
                atm.vrt_i = vrt_i
                atm.root = i
                atm.ring_index = rnum
                
                if (gp_i == i): # near-root corr
                    atm.groot = (par_i,nodes[2])
                else:
                    atm.groot = (par_i,gp_i)
                mol.vatms.append(atm)

def is_biaryl_ring(mol,ring1,ring2):
    '''Check conjugation between two rings'''
    # first check if shares any atom
    for atm1 in ring1.atms:
        for atm2 in ring2.atms:
            if atm1 == atm2:
                return False
    # Then figure out if any atm pair connected
    for atm1 in ring1.atms:
        for atm2 in ring2.atms:
            for bond in mol.bonds:
                if (atm1,atm2) != (bond.atm1,bond.atm2) and \
                   (atm1,atm2) != (bond.atm2,bond.atm1):
                    continue

                #check if atm1,atm2 are connected by ring
                is_connected_by_ring = False
                for ring in mol.rings:
                    if atm1 in ring.atms and atm2 in ring.atms:
                        is_connected_by_ring = True
                        break
                if not is_connected_by_ring:
                    return (atm1,atm2)
    return False

def define_conjugation(mol):
    '''Defines conjugations in bonds in Molecule based on heuristics'''
    # 1. Rotable ring-ring
    mol.biaryl_rings = []
    mol.biaryl_pivots = []
    if len(mol.rings) > 1:
        for i,ring1 in enumerate(mol.rings[:-1]):
            for j,ring2 in enumerate(mol.rings[i+1:]):
                biaryl_pivot = is_biaryl_ring(mol,ring1,ring2)
                if biaryl_pivot:
                    mol.biaryl_rings.append([i,i+j+1])
                    mol.biaryl_pivots.append(biaryl_pivot)
        if mol.option.verbose:
            print( 'Added regular biaryl-pivots: ')
            for a1,a2 in mol.biaryl_pivots:
                print( " (%s,%s)"%(mol.atms[a1].name,mol.atms[a2].name))

    ## 2. Pseudo ring-ring: non-ring connected to ring
    if len(mol.rings) > 0:
        biaryl_pivot_extra = []
        for i,ring1 in enumerate(mol.rings_aro):
            biaryl_pivot_extra += search_special_biaryl_ring(mol,ring1)
        if mol.option.verbose:
            print( 'Added extra biaryl-pivots: ')
            for a1,a2 in biaryl_pivot_extra:
                print( " (%s,%s)"%(mol.atms[a1].name,mol.atms[a2].name))
        mol.biaryl_pivots_extra = biaryl_pivot_extra
        mol.biaryl_pivots += biaryl_pivot_extra

    # 3. Amide bond; trick to allow rotation around amide bonds
    for i,bond in enumerate(mol.bonds):
        aclass1 = ACLASS_ID[mol.atms[bond.atm1].aclass]
        aclass2 = ACLASS_ID[mol.atms[bond.atm2].aclass]
        is_amide_connection = (aclass1 in ['Nad','Nad3']) and (aclass2 in ['Nad','Nad3'])
        if is_amide_connection: # append into special biaryl
            mol.biaryl_pivots.append((bond.atm1,bond.atm2))
            
    if mol.option.opt.reassign_biaryl_aclass:
        reassign_biaryl_atypes(mol)
    else:
        assign_bond_conjugation(mol)

    if mol.option.verbose:
        print('Atom_Puckering: ', [mol.atms[atm].name for atm in mol.atms_puckering])
        print('Atom_Aro: ', [mol.atms[atm].name for atm in mol.atms_aro])
        print('Rings: ', [[mol.atms[atm].name for atm in ring.atms] for ring in mol.rings])
        if len(mol.biaryl_rings) > 0 :
            print('BiarylRings:', mol.biaryl_rings, ' BiarylAxes: ') 
            for a1,a2 in mol.biaryl_pivots:
                print((mol.atms[a1].name,mol.atms[a2].name))
            print()

# Deprecated: For atom-type-based torsion assignment
def reassign_biaryl_atypes(mol):
    # reassign atype
    # Turned off as of Oct 2018 with new bond-based torsion assignment logic
    for a1,a2 in mol.biaryl_pivots:
        if mol.atms[a1].aclass in [ACLASS_ID.index('CR'),ACLASS_ID.index('CRp'),
                                   ACLASS_ID.index('CD1'),ACLASS_ID.index('CD2'),ACLASS_ID.index('CD'),ACLASS_ID.index('CDp')]: #carbon
            mol.atms[a1].aclass = ACLASS_ID.index('CRb')
        elif mol.atms[a1].atype == ATYPES.index('N'):
            mol.atms[a1].aclass = ACLASS_ID.index('NGb')
            
        if mol.atms[a2].aclass in [ACLASS_ID.index('CR'),ACLASS_ID.index('CRp'),
                                   ACLASS_ID.index('CD1'),ACLASS_ID.index('CD2'),ACLASS_ID.index('CDp'),ACLASS_ID.index('CDp')]: #carbon
            mol.atms[a2].aclass = ACLASS_ID.index('CRb')
        elif mol.atms[a2].atype == ATYPES.index('N'):
            mol.atms[a2].aclass = ACLASS_ID.index('NGb')

# For bond-type-based torsion assignment
# use bond-order + biaryl_pivot info
def assign_bond_conjugation(mol):    
    for ibond,bond in enumerate(mol.bonds):
        # Consider not conjugated if any of two is sp3
        if mol.atms[bond.atm1].hyb == 3 or mol.atms[bond.atm2].hyb == 3: continue

        # Jan 2019: To avoid forcing planar forms for puckering rings
        # Also not conjugated if belongs to a puckering ring
        bond_at_puckering_ring = False
        for ring in mol.rings:
            if ring.type < 4 and ring.has((bond.atm1,bond.atm2)):
                bond_at_puckering_ring = True
                break
        if bond_at_puckering_ring: continue

        # special conjugation through ring(N=C)-N-H
        # THIS PART SHOULD BE REVISITED WITH "RING" BONDTYPE IN FUTURE
        is_ring_NCNH = False
        if (bond.atm1 in mol.atms_aro and bond.atm2 not in mol.atms_aro) or \
           (bond.atm2 in mol.atms_aro and bond.atm1 not in mol.atms_aro):
            # i: ring-member connected to other branch
            # j: non-ring-member connected to i
            if (bond.atm1 in mol.atms_aro and bond.atm2 not in mol.atms_aro):
                i,j = (bond.atm1,bond.atm2)
            else:
                i,j = (bond.atm2,bond.atm1)
                
            type_i = ATYPES[mol.atms[i].atype]
            type_j = ATYPES[mol.atms[j].atype]
            if mol.atms[j].has_H and (type_i,type_j) == ('C','N'):
                for k,dumm in mol.atms[i].bonds:
                    if mol.in_same_ring(i,k) == 0 or k == j: continue
                    if (ATYPES[mol.atms[k].atype] == 'N') and (mol.atms[k].hyb != 3):
                        is_ring_NCNH = True
                        break

        if is_ring_NCNH:
            print("Torsion around %s-%s assigned as conjugated by [ring N=C]-N-H rule"%(mol.atms[bond.atm1].name,mol.atms[bond.atm2].name))
            mol.bonds[ibond].is_conjugated = True
            continue

        # also skip if is part of >= 7-membered ring (whether aromatic or not)
        # because conjugation start to break (more influenced by neighbors than conjugation)
        # THIS PART SHOULD BE REVISITED WITH "RING" BONDTYPE IN FUTURE
        minimum_ring_size = 100
        for ring in mol.rings:
            if ring.has((bond.atm1,bond.atm2)):
                if ring.natms < minimum_ring_size:
                    minimum_ring_size = ring.natms
        if minimum_ring_size < 100 and minimum_ring_size > 6:
            continue

        # simpler logic over-predicts conjugated rings
        #if mol.in_same_ring(bond.atm1,bond.atm2) == 1: continue

        #print( "assign conjugation: ", mol.atms[bond.atm1].name, mol.atms[bond.atm2].name,
        #       (bond.atm1,bond.atm2) in mol.biaryl_pivots or \
        #       (bond.atm2,bond.atm1) in mol.biaryl_pivots,
        #       ACLASS_ID[mol.atms[bond.atm1].aclass] in CONJUGATING_ACLASSES,
        #       ACLASS_ID[mol.atms[bond.atm2].aclass] in CONJUGATING_ACLASSES)

        if ACLASS_ID[mol.atms[bond.atm1].aclass] not in CONJUGATING_ACLASSES or \
           ACLASS_ID[mol.atms[bond.atm2].aclass] not in CONJUGATING_ACLASSES: continue

        if ((bond.atm1,bond.atm2) not in mol.biaryl_pivots) and \
           ((bond.atm2,bond.atm1) not in mol.biaryl_pivots):
            mol.bonds[ibond].is_conjugated = True

# Define "CHI"s
def define_rotable_torsions(mol):
    mol.chiatms = []
    mol.chitypes = []
    mol.chiextra = ""

    hapol_torsion_id = []
    hpol_torsion_type = [0 for k in mol.torsions]

    # 1. First get list of hydrogen torsions
    for i,torsion in enumerate(mol.torsions):
        aclasses = [ACLASS_ID[mol.atms[atm].aclass] for atm in torsion]
        # Apolar hydrogen torsions
        if (aclasses[0] in ACLASS_HAPOL) or (aclasses[3] in ACLASS_HAPOL):
            hapol_torsion_id.append(i)

        # PolarH torsions
        elif (aclasses[0] in ACLASS_HPOL) and (aclasses[3] not in ACLASS_HPOL):
            stem = mol.atms[torsion[1]]
            if stem.hyb == 2:
                if len(stem.bonds) < 3:
                    hpol_torsion_type[i] = 2
            elif stem.hyb == 3: #mol.torsion_orders[i] == 1:
                hpol_torsion_type[i] = 3

        elif (aclasses[3] in ACLASS_HPOL) and (aclasses[0] not in ACLASS_HPOL):
            stem = mol.atms[torsion[2]]
            if stem.hyb == 2:
                if len(stem.bonds) < 3:
                    hpol_torsion_type[i] = 2
            elif stem.hyb == 3: #mol.torsion_orders[i] == 1:
                hpol_torsion_type[i] = 3

    # 1-2. Get estimations of num_H_conf
    num_H_confs = 1
    covered = []
    for i,atms in enumerate(mol.torsions):
        if (atms[1],atms[2]) in covered or (atms[2],atms[1]) in covered: continue
        if hpol_torsion_type[i] == 2: num_H_confs *= 6
        elif hpol_torsion_type[i] == 3: num_H_confs *= 9
        covered.append((atms[1],atms[2]))
        
    #print 'Total num_H_conf: ', num_H_confs
    if num_H_confs <= mol.max_confs: mol.chiextra = "1 20"
    else: mol.chiextra = "0"

    # 2. CHI assignment; read-in from mol.option
    # mol.torsions SHOULD have been ordered such that
    # heavyatom-only comes first followed by H-containing ones
    ring_cuts = []
    for ring in mol.rings: ring_cuts.append(ring.cut_bond)

    covered = []
    for i,atms in enumerate(mol.torsions):
        atm0 = mol.atms[atms[0]]
        atm1 = mol.atms[atms[1]]
        atm2 = mol.atms[atms[2]]
        atm3 = mol.atms[atms[3]]
        border = mol.bond_order(atms[1],atms[2])

        if (atms[1],atms[2]) in covered or (atms[2],atms[1]) in covered: continue #avoid
        if ((atm0.root != atms[1]) and (atm1.root != atms[0])) or ((atm2.root != atms[3]) and (atm3.root != atms[2])): continue
            
        if atm1.root == atms[2]:
            atms_ordered = [atms[3],atms[2],atms[1],atms[0]] #last is tipatm
        elif atm2.root == atms[1]:
            atms_ordered = [atms[0],atms[1],atms[2],atms[3]] #last is tipatm
        else:
            continue

        # avoid if any of bonds defined as cut_bond
        if (atms[0],atms[1]) in ring_cuts or (atms[1],atms[0]) in ring_cuts or \
           (atms[1],atms[2]) in ring_cuts or (atms[2],atms[1]) in ring_cuts or \
           (atms[2],atms[3]) in ring_cuts or (atms[3],atms[2]) in ring_cuts:
            continue
        
        # why did I put this logic here...?
        # check if the torsion is non-ATorder but connected to ATorder
        FT_connected = False
        for j in range(len(atms)-1):
            a1,a2 = atms[j],atms[j+1]
            if a1 in mol.ATorder and a2 in mol.ATorder: continue
            if (mol.atms[a1].root not in [a2,a1]) and (mol.atms[a2].root not in [a1,a2]):
                FT_connected = True
                break
        if FT_connected: continue

        covered.append((atms[1],atms[2]))
    
        is_biaryl_pivot = ((atms[1],atms[2]) in mol.biaryl_pivots) or \
                          ((atms[2],atms[1]) in mol.biaryl_pivots)

        # Bug fix with biaryl not defined as CHI: Oct17 2018 
        if ((atms[1] in mol.atms_aro) and (atms[2] in mol.atms_aro)):
            if border>1: continue
            same_ring_order = mol.in_same_ring(atms[1],atms[2])
            if same_ring_order == 4:
                continue # also skip if belongs to the same "aromatic" ring
            elif same_ring_order == 0: #not in same ring
                if not mol.option.opt.report_ringring_chi and is_biaryl_pivot:
                    continue

        # below are optional
        # Apolar hydrogen chis
        if (not mol.option.opt.report_Hapol_chi) and (i in hapol_torsion_id): continue 
    
        if (mol.bond_order(atms[1],atms[2]) > 1):
            atype1 = ACLASS_ID[atm1.aclass]
            atype2 = ACLASS_ID[atm2.aclass]
            is_amide_bond = ((atype1 == 'Nad') and (atype2 in ['CDp','CRb'])) or \
                            ((atype2 == 'Nad') and (atype1 in ['CDp','CRb']))

            #is_aliphatic_bond = ((atype1 == 'CD1') and (atype2 == 'CD1') and \
            #                     (atms[1] not in mol.atms_ring) and (atms[2] not in mol.atms_ring))

            #print atype1, atype2, is_amide_bond
            #print( mol.atms[atms[1]].name,mol.atms[atms[2]].name,is_aliphatic_bond)
            #if is_aliphatic_bond:
            #    pass
            if (mol.option.opt.report_amide_chi):
                if (not is_amide_bond): continue
            elif (not mol.option.opt.report_nbonded_chi and not is_biaryl_pivot):
                continue

        # skip conjugated polarH chi
        elif mol.bond_conjugated(atms[1],atms[2]):
            nH1 = 0
            for j,order in atm1.bonds:
                if mol.atms[j].is_H: nH1+=1
            nH2 = 0
            for j,order in atm2.bonds:
                if mol.atms[j].is_H: nH2+=1

            if nH1 == len(atm1.bonds)-1 or nH2 == len(atm2.bonds)-1:
                continue

        if (not mol.option.opt.report_puckering_chi) and \
           (atms[1] in mol.atms_puckering) and (atms[2] in mol.atms_puckering):
            continue

        mol.chiatms.append(atms_ordered)
        
        if hpol_torsion_type[i] == 2: #sp2
            mol.chitypes.append('sp2')
        elif hpol_torsion_type[i] == 3: #sp3
            mol.chitypes.append('sp3')
        elif i in hapol_torsion_id: #sp3
            mol.chitypes.append('sp3H')
        #elif (atms[1] in mol.atms_puckering) and (atms[2] in mol.atms_puckering):
        #    mol.chitypes.append('pucker') #turn this functionality for now... add separate flag for this behavior
        else:
            mol.chitypes.append('')
