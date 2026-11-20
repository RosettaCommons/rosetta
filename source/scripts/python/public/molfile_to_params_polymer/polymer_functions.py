#!/usr/bin/env python

import os, sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

# Magic spell to make sure Rosetta python libs are on the PYTHONPATH:
#sys.path.append( os.path.abspath( sys.path[0] ) )
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

#from rosetta_py.io.mdl_molfile import *
#from rosetta_py.utility.rankorder import argmin
from rosetta_py.utility import r3
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

import functools

def polymer_assign_backbone_atom_types(m):
    # first get POLY flags from molfile
    for line in m.footer:
        if line.startswith("M  POLY_N_BB"):
            num = atom_num( line.split()[2], m )
            m.atoms[ num ].poly_backbone = True
            m.atoms[ num ].poly_n_bb = True
        if line.startswith("M  POLY_CA_BB"):
            num = atom_num( line.split()[2], m )
            m.atoms[ num ].poly_backbone = True
            m.atoms[ num ].poly_ca_bb = True
        if line.startswith("M  POLY_O_BB"):
            num = atom_num( line.split()[2], m )
            m.atoms[ num ].poly_backbone = True
            m.atoms[ num ].poly_o_bb = True
        if line.startswith("M  POLY_C_BB"):
            num = atom_num( line.split()[2], m )
            m.atoms[ num ].poly_backbone = True
            m.atoms[ num ].poly_c_bb = True
        if line.startswith("M  POLY_UPPER"):
            num = atom_num( line.split()[2], m )
            m.atoms[ num ].poly_backbone = True
            m.atoms[ num ].poly_upper = True
        if line.startswith("M  POLY_LOWER"):
            num = atom_num( line.split()[2], m )
            m.atoms[ num ].poly_backbone = True
            m.atoms[ num ].poly_lower = True
        if line.startswith("M  POLY_BACKBONE"):
            for num in [ atom_num(s, m) for s in line.split()[2:] ]:
                m.atoms[ num ].poly_backbone = True

def polymer_assign_backbone_atom_names(atoms, bonds, peptoid):
    ''' modifies the rosetta atom type of the hydorgen attached to the alpha carbon and backbone nitrogen to be the correct
    rosetta atom types '''
    # heavy atom atom names
    for atom in atoms:
        if atom.poly_ca_bb:
            atom.ros_type = "CAbb"
            atom.pdb_name = " CA "
        elif atom.poly_n_bb:
            if atom.ros_type == "":
                atom.ros_type = "Nbb"
            atom.pdb_name = " N  "
            if peptoid:
                atom.mm_type = "NXX"
        elif atom.poly_o_bb:
            atom.ros_type = "OCbb"
            atom.pdb_name = " O  "
        elif atom.poly_c_bb:
            atom.ros_type = "CObb"
            atom.pdb_name = " C  "
        # Upper & Lower handled in a different function.
    # alpha carbon hydrogen(s)
    ca_h_bonds = []
    index = 0;
    for bond in bonds:
        if bond.a1.poly_ca_bb and bond.a2.is_H :
            bond.a2.ros_type = "Hapo"
            ca_h_bonds.append( index)
            if not peptoid:
                bond.a2.pdb_name = " HA "
        elif bond.a1.is_H and bond.a2.poly_ca_bb :
            bond.a1.ros_type = "Hapo"
            ca_h_bonds.append( index)
            if not peptoid:
                bond.a1.pdb_name = " HA "
        index = index + 1
    # for the special case of Glycine or peptoid
    if len(ca_h_bonds) > 1:
        bond_id = 0
        for index in ca_h_bonds :
            if bonds[ index].a1.is_H:
                bonds[ index].a1.pdb_name = "%dHA " % (bond_id + 1)
            else :
                bonds[ index].a2.pdb_name = "%dHA " % (bond_id + 1)
            bond_id = bond_id + 1

    # backbone nitrogen hydrogen(s)
    for bond in bonds:
        if bond.a1.poly_n_bb and bond.a2.is_H :
            bond.a2.ros_type = "HNbb"
            bond.a2.pdb_name = " H  "
        elif bond.a1.is_H and bond.a2.poly_n_bb :
            bond.a1.ros_type = "HNbb"
            bond.a1.pdb_name = " H  "

def polymer_assign_connection_atom_names(atoms, bonds, peptoid):
    ''' modifies the names of the upper and lower connections '''
    for atom in atoms:
        if atom.poly_upper:
            atom.ros_type = "X"
            atom.pdb_name = "UPPER"
        elif atom.poly_lower:
            atom.ros_type = "X"
            atom.pdb_name = "LOWER"


def polymer_assign_ignored_atoms_bonds(m):
    ''' sets the ignore boolean for each atom in the list and for each bond with at least one atom in the list'''
    ignore_list = []
    for line in m.footer:
        if line.startswith("M  POLY_IGNORE"):
            ignore_list = [ atom_num(d, m) for d in line.split()[2:] ]
    #atoms
    for i,a in enumerate(m.atoms):
        if i in ignore_list:
            a.poly_ignore = True
    # bonds
    print("the molecule has %d atoms" % len(m.atoms))
    for bond in m.bonds:
        for i in ignore_list:
            #print("ignore atom: %d" % i)
            if m.atoms[i] == bond.a1 or m.atoms[i] == bond.a2:
                bond.poly_ignore = True

def polymer_assign_pdb_like_atom_names_to_sidechain(atoms, bonds, peptoid):
    ''' Assign PDB like names to atoms. PDB names are based on the path distance from the alpha carbon
    Greek Alphabet (it has been extended to lower-case letters for large non-canonical AAs like lanthanide-binding tags): Alpha, Beta, Gamma, Delta, Epsilon, Zeta, Eta, Theta, Iota, Kappa, Lambda, Mu, Nu, Xi, Omicron, Pi, Rho, Sigma, Tau, Upsilon, Phi, Chi, Psi, Omega, alpha, beta, gamma, delta, epsilon, zeta, eta, theta, iota, kappa, lambda, mu, nu, xi, omicron, pi, rho, sigma, tau, upsilon, phi, chi, psi, omega'''
    # greek alphabet eta, tao and omega are skipped because they are the same as previous letters
    greek_alphabet = ['A', 'B', 'G', 'D', 'E', 'Z', 'T', 'I', 'K', 'L', 'M', 'N', 'X', 'O', 'P', 'R', 'S', 'U', 'P', 'C', 'a', 'b', 'g', 'd', 'e', 'z', 't', 'i', 'k', 'l', 'm', 'n', 'x', 'o', 'p', 'r', 's', 'u', 'p', 'c']
    def get_atom_num(elem):
        elem_atom_num = {'B':5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'NA': 11, 'MG': 12, 'P':15, 'S':16, 'CL':17, 'K':19, 'CA':20,
                         'FE':26, 'ZN':30, 'BR':35, 'I':53, 'SE':34, 'X':-1}
        try:
            return elem_atom_num[elem]
        except:
            print("WARNING: element type %s not recognized, pdb naming may be incorrect"%elem)
            return 99

    # find alpha carbon or the ""alpha nitrogen" for peptoids (still called ca_index below)
    #print "PEPTOID" , peptoid
    if peptoid:
        for ca_index, atom in enumerate(atoms):
            if atom.poly_n_bb:
                #print atom
                break
    else:
        for ca_index, atom in enumerate(atoms):
            if atom.poly_ca_bb:
                break
    #print "CA or NA index is %d" % ca_index
    # assign heavy atom and hydrogen pdb_elem
    for atom in atoms:
        atom.pdb_elem = atom.elem
    # assign heavy atom pdb_greek_dist
    def path_dist(a,b):
        return 1
    na = len(atoms) # Number of Atoms
    all_all_dist = [ [1e100] * na for i in range(na) ]
    nbrs = dict([ (a,set()) for a in atoms ])
    for a in atoms:
        #nbrs[a].update([b.a2 for b in a.bonds if b.a2.is_H == False and b.a2.poly_ignore == False and b.a2.poly_backbone == False])
        nbrs[a].update([b.a2 for b in a.bonds if b.a2.is_H == False and b.a2.poly_ignore == False and b.a2.poly_n_bb == False and b.a2.poly_c_bb == False and b.a2.poly_o_bb == False and b.a2.poly_upper == False and b.a2.poly_lower == False])
    for i in range(0,na):
        all_all_dist[i] = dijkstra( start = atoms[i], nodes = atoms, nbr = lambda a: nbrs[a], dist = path_dist )
    #print "ALL TO ALL DIST"
    #debug
    #for i,a in enumerate(atoms):
        #print a, all_all_dist[i]
    #print "ALL TO ALL DIST CA INDEX"
    #print all_all_dist[ca_index] #DEBUG
    for i, a in enumerate(atoms):
        if peptoid:
            if not a.is_H and not a.poly_ignore and not a.poly_n_bb and not a.poly_c_bb and not a.poly_o_bb and not a.poly_upper and not a.poly_lower:
                #print "ATOM: ", a
                #print "DISTANCE: %f" % (all_all_dist[ca_index][i]-1)
                a.pdb_greek_dist = greek_alphabet[all_all_dist[ca_index][i]-1]
        else:
            if not a.is_H and not a.poly_ignore and not a.poly_backbone:
                #print "ATOM: ", a
                #print "DISTANCE: %d" % all_all_dist[ca_index][i]
                a.pdb_greek_dist = greek_alphabet[all_all_dist[ca_index][i]]
    #debug = [a.pdb_greek_dist for a in atoms if not a.is_H ] #DEBUG
    #print debug #DEBUG
    # assign heavy atom pdb_postfix_num (stupidly inefficient)
    def compare_atom_num(x,y):
        if get_atom_num(atoms[x].elem) > get_atom_num(atoms[y].elem):
            return -1
        elif get_atom_num(atoms[x].elem) == get_atom_num(atoms[y].elem):
            heavy_x = [b for b in atoms[x].bonds if not b.a2.is_H]
            heavy_y = [b for b in atoms[y].bonds if not b.a2.is_H]
            if len(heavy_x) > len(heavy_y):
                return -1
            elif len(heavy_x) == len(heavy_y):
                return 0
            elif len(heavy_x) < len(heavy_y):
                return 1
        elif get_atom_num(atoms[x].elem) < get_atom_num(atoms[y].elem):
            return 1
    for i, g in enumerate(greek_alphabet):
        temp = [j for j,a in enumerate(atoms) if a.pdb_greek_dist == greek_alphabet[i]]
        if len(temp) > 1 or (len(temp) == 1 and len(atoms[temp[0]].bonds) == 1):
            temp.sort(key=functools.cmp_to_key(compare_atom_num))
            used_index = []
            #Looping twice, once to assign easy parent-based postfixes, once to assign all others
            for t in temp:
                attached_num = [b.a2.pdb_postfix_num for b in atoms[t].bonds if b.a2.pdb_postfix_num != " " and greek_alphabet.index(b.a2.pdb_greek_dist) < i]
                #If no numbers come, different numbers come in, or the number is already taken, assign a number normally later
                if len(attached_num) == 1 and attached_num[0] not in used_index:
                    atoms[t].pdb_postfix_num = attached_num[0]
                    used_index.append(attached_num[0])
            k = 1
            for t in temp:
                while str(k) in used_index:
                    k += 1
                if atoms[t].pdb_postfix_num == " ":
                    atoms[t].pdb_postfix_num = "%d" % k
                    used_index.append(str(k))

        #CAA's label bonded atoms with the same postfix number, even if there's only one
        elif len(temp) == 1:
            t = temp[0]
            attached_num = [b.a2.pdb_postfix_num for b in atoms[t].bonds if b.a2.pdb_postfix_num != " " and greek_alphabet.index(b.a2.pdb_greek_dist) < i]
            #If you have different numbers coming in, don't number it at all
            if len(attached_num) == 1:
                atoms[t].pdb_postfix_num = attached_num[0]

    #debug
    #for a in atoms:
        #if a.poly_ca_bb:
            #print "DEBUG CA_BB: ", a, ":", a.pdb_prefix_num, ":", a.pdb_elem, ":", a.pdb_greek_dist, ":", a.pdb_postfix_num
            #for b in a.bonds:
                #print b
    # assign hydrogen pdb_greek_dist and pdb_postfix_num
    for a in atoms:
        if a.is_H and not a.poly_c_bb:
            a.pdb_greek_dist = a.bonds[0].a2.pdb_greek_dist
            a.pdb_postfix_num = a.bonds[0].a2.pdb_postfix_num
    # assign hydrogen pdb_prefix_num
    if peptoid:
        for a in atoms:
            if not a.is_H and not a.poly_ignore and not a.poly_n_bb and not a.poly_c_bb and not a.poly_o_bb and not a.poly_upper and not a.poly_lower:
                attached_h = [atoms.index(b.a2) for b in a.bonds if b.a2.is_H == True]
                if len(attached_h) > 1:
                    for i,ah in enumerate(attached_h):
                        blah = i + 1
                        atoms[ah].pdb_prefix_num = "%d" % blah
    else:
        for a in atoms:
            if not a.is_H and not a.poly_backbone and not a.poly_ignore:
                attached_h = [atoms.index(b.a2) for b in a.bonds if b.a2.is_H == True]
                if len(attached_h) > 1:
                    for i,ah in enumerate(attached_h):
                        blah = i + 1
                        atoms[ah].pdb_prefix_num = "%d" % blah
    # assign full pdb name
    for a in atoms:
        a.pdb_name = a.pdb_prefix_num + a.pdb_elem + a.pdb_greek_dist + a.pdb_postfix_num
    #debug
    #for a in atoms:
        #if a.poly_ca_bb:
            #for b in a.bonds:
                #print "DEBUG: ", b.a2.pdb_name

def polymer_reorder_atoms(molfile):
    ''' Reorders the atoms acording to the pdb ordering so that the order of the internal coords is correct '''
    def poly_atom_cmp(atom1, atom2):
        ''' Sorts based on special polymer backbone type, greek letter distance, postfix num, prefix num '''
        greek_alphabet = { ' ':0, 'A':1, 'B':2, 'G':3, 'D':4, 'E':5, 'Z':6, 'T':7, 'I':8, 'K':9, 'L':10, 'M':11, 'N':12, 'X':13, 'O':14, 'P':15, 'R':16, 'S':17, 'U':18, 'P':19, 'C':20, 'W':21, 'W':22, 'W':23, 'W':24, 'W':25, 'W':26, 'W':27, 'W':28, 'W':29, 'W':30, 'W':31, 'a':32, 'b':33, 'g':34, 'd':35, 'e':36, 'z':37, 't':38, 'i':39, 'k':40, 'l':41, 'm':42, 'n':43, 'x':44, 'o':45, 'p':46, 'r':47, 's':48, 'u':49, 'p':50, 'c':51 }
        # ignore
        if atom1.poly_ignore == True and atom2.poly_ignore == True: return 0
        elif atom1.poly_ignore == True and atom2.poly_ignore == False: return 1
        elif atom1.poly_ignore == False and atom2.poly_ignore == True: return -1
        elif atom1.poly_ignore == False and atom2.poly_ignore == False:
                    # special poly types
                    if atom1.poly_backbone == True and atom2.poly_backbone == False: return -1
                    elif atom1.poly_backbone == False and atom2.poly_backbone == True: return 1
                    elif atom1.poly_backbone == True and atom2.poly_backbone == True:
                        if atom1.poly_n_bb == True: return -1
                        elif atom1.poly_ca_bb == True and atom2.poly_n_bb != True: return -1
                        elif atom1.poly_c_bb == True and atom2.poly_n_bb != True and atom2.poly_ca_bb != True: return -1
                        elif atom1.poly_o_bb == True and atom2.poly_n_bb != True and atom2.poly_ca_bb != True and atom2.poly_c_bb != True: return -1
                        else: return 1
                    elif atom1.poly_backbone == False and atom2.poly_backbone == False:
                        # upper
                        if atom1.poly_upper == True and atom2.poly_upper == True: return 0 # this should never happen
                        if atom1.poly_upper == True and atom2.poly_upper == False: return 1
                        elif atom1.poly_upper == False and atom2.poly_upper == True: return -1
                        elif atom1.poly_upper == False and atom2.poly_upper == False:
                            # lower
                            if atom1.poly_lower == True and atom2.poly_lower == True: return 0 # this should never happen
                            if atom1.poly_lower == True and atom2.poly_lower == False: return 1
                            elif atom1.poly_lower == False and atom2.poly_lower == True: return -1
                            elif atom1.poly_lower == False and atom2.poly_lower == False:
                                # hydrogen
                                if atom1.elem == 'H' and atom2.elem != 'H': return 1
                                # greek distance
                                if greek_alphabet[atom1.pdb_greek_dist] < greek_alphabet[atom2.pdb_greek_dist]: return -1
                                elif greek_alphabet[atom1.pdb_greek_dist] > greek_alphabet[atom2.pdb_greek_dist]: return 1
                                elif greek_alphabet[atom1.pdb_greek_dist] == greek_alphabet[atom2.pdb_greek_dist]:
                                    # postfix num
                                    if atom1.pdb_postfix_num < atom2.pdb_postfix_num: return -1
                                    elif atom1.pdb_postfix_num > atom2.pdb_postfix_num: return 1
                                    elif atom1.pdb_postfix_num == atom2.pdb_postfix_num:
                                        # prefix num
                                        if atom1.pdb_prefix_num < atom2.pdb_prefix_num: return -1
                                        elif atom1.pdb_prefix_num > atom2.pdb_prefix_num: return 1
                                        else: return 0 #
    molfile.atoms.sort(key=functools.cmp_to_key(poly_atom_cmp) )
