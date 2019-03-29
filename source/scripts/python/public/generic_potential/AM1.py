#!/usr/bin/env python3
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
'''
Python implementation of AM1 partial charge calculation.

Author: Frank DiMaio 
'''
import os, sys

import math
import json
import io

import numpy
import numpy.linalg
import numpy.testing
import scipy.sparse.csgraph

from numba import autojit, jit, prange

import AM1_overlap as ovr

direc = os.path.dirname( os.path.abspath(__file__) )
#MYFILE = os.path.abspath(__file__)
#direc = MYFILE.replace('AM1.py','')

#from concurrent import futures
# multithreading params
#nthreads = 4

# constants
BohrRadius = 0.529177  # in Angstroms
Hartree = 27.21138  # in eV
EVtoKcal = 23.061

# database
#print direc, "%s/am1params.json"%direc
AM1db = json.load(open("%s/am1params.json"%direc, "r"))

# per-orbital data
orbital_dtype = numpy.dtype([
    ('atomno', numpy.int),
    ('type', numpy.int),   # ['s','px','py','pz']
    ('U' , numpy.float),
    ('alpha' , numpy.float),
    ('beta' , numpy.float),
])

# per-atom data
atom_dtype = numpy.dtype([
    ('name' , numpy.str, 4),
    ('elt' , numpy.str, 2),
    ('X' , numpy.float, 3),
    ('B' , numpy.float, 1),
    ('atno' , numpy.int),  # atomic number
    ('zetaS' , numpy.float),
    ('zetaP' , numpy.float)
])

# epsilon data
epsilon_dtype = numpy.dtype([
    ('ssss' , numpy.str, 4),
    ('ssspx' , numpy.str, 2),
    ('ssspy' , numpy.float, 3),
])


@jit(nopython=True,fastmath=True)
def cross(a, b):
    result = numpy.zeros(3)
    result[0] = a[1] * b[2] - a[2] * b[1]
    result[1] = a[2] * b[0] - a[0] * b[2]
    result[2] = a[0] * b[1] - a[1] * b[0]
    return result

### jitted get Rs: atom-pair distances
@jit(nopython=True,fastmath=True)
def inline_Rs(iatm, jatm, atoms, Rvecs, Rs, frame):
    def tri_idx(i,j,N):
        return int( (N-1-(i+1)/2)*i + j - 0.5)

    natoms = atoms.shape[0]
    for i in range(len(iatm)):
        a_i = iatm[i]
        a_j = jatm[i]
        if (a_i > a_j):
            a_j = iatm[i]
            a_i = jatm[i]

        Rvecs[a_i,a_j,:] = atoms[a_i,:]-atoms[a_j,:]
        Rs[a_i,a_j] = numpy.sqrt( numpy.dot( Rvecs[a_i,a_j,:], Rvecs[a_i,a_j,:] ))
        Rvecs[a_i,a_j,:] /= Rs[a_i,a_j]
        Rs[a_i,a_j] /= BohrRadius

        Rvecs[a_j,a_i,:] = -Rvecs[a_i,a_j,:]
        Rs[a_j,a_i] = Rs[a_i,a_j]

        # frame
        Z = Rvecs[a_j,a_i,:]
        X = numpy.copy( Z )
        if (abs(X[2]) < 0.7):
            X[0],X[1],X[2] = Z[1],-Z[0],Z[2]
        else:
            X[0],X[1],X[2] = Z[2],-Z[1],Z[0]
        X -= numpy.sum(X*Z)*Z
        X /= numpy.sqrt( numpy.dot( X, X ))
        Y = cross(Z,X)

        idx1 = tri_idx(a_i,a_j,natoms)
        frame[idx1,0,0]=1
        frame[idx1,1,1:]=X
        frame[idx1,2,1:]=Y
        frame[idx1,3,1:]=Z

### jitted get V: the one-electron matrix
@jit(nopython=True,fastmath=True)
def inline_V(V, Zs, epsilon, orbatm, orbnum, norbs, idxorbs):
    def tri_idx(i,j,N):
        return int( (N-1-(i+1)/2)*i + j - 0.5)

    norbstot = V.shape[0]
    natms = idxorbs.shape[0]

    for atm_i in range(natms):
        for orb_i1 in range(idxorbs[atm_i],idxorbs[atm_i]+norbs[atm_i]):
            for orb_i2 in range(orb_i1,idxorbs[atm_i]+norbs[atm_i]):
                Vi = 0
                for atm_j in range(atm_i):
                    idx_ij = tri_idx(atm_j,atm_i,natms)
                    Vi -= Zs[atm_j] * epsilon[idx_ij, orbnum[orb_i1], orbnum[orb_i2], 0, 0]

                for atm_j in range(atm_i+1,natms):
                    idx_ij = tri_idx(atm_i,atm_j,natms)
                    Vi -= Zs[atm_j] * epsilon[idx_ij, 0, 0, orbnum[orb_i1], orbnum[orb_i2]]

                V[orb_i1,orb_i2] = Vi
                V[orb_i2,orb_i1] = Vi

### jitted get_epsilon: two electron two center integrals
@jit(nopython=True,parallel=True,fastmath=True)
def inline_epsilon(atompairs, P1s, P2s, orb_is, orb_js, atmj, atmi, Rs, Rvecs, frame):
    def f1(r, rhoA, rhoB):
        return 1.0 / numpy.sqrt(
            numpy.square(r) + numpy.square(rhoA + rhoB))

    def dist(X, Y, Rz):
        r = numpy.sqrt(
            numpy.square(Y[0] - X[0]) +
            numpy.square(Y[1] - X[1]) +
            numpy.square(Y[2] - X[2] + Rz))
        return r

    def mp_mp(P1, P2, Rz):
        J = 0
        for i in range(4):
            for j in range(4):
                J += P1[i,4] * P2[j,4] * f1(dist(P1[i,:3], P2[j,:3], Rz), P1[i,3], P2[j,3])
        return J * Hartree

    natmpairs = atompairs.shape[0]
    epsilon = numpy.zeros( (natmpairs, 4,4,4,4), dtype=numpy.float64 )

    # calc in canonic frame
    for iidx in prange(natmpairs):
        i = atompairs[iidx]
        scratch = numpy.zeros( (10, 10), dtype=numpy.float64 )
        #scratch.fill(0)
        epsilon_i = epsilon[iidx,...]

        for j in range(10):
            for k in range(10):
                J_ijk = mp_mp( P1s[i,j,:,:], P2s[i,k,:,:], Rs[atmj[i], atmi[i]] )

                if abs(J_ijk) > 1e-10:
                    scratch[j,k] = J_ijk

        #idxI,idxJ = numpy.where(scratch != 0.0)

        for j1 in range(4):
            for j2 in range(j1,4):
                for k1 in range(4):
                    for k2 in range(k1,4):
                        J = 0

                        #for N in range(len(idxI)):
                        for ind_i in range(10):
                            for ind_j in range(10):
                                if (scratch[ind_i,ind_j] == 0):
                                    continue

                                orb_i1 = orb_is[ind_i]
                                orb_i2 = orb_js[ind_i]
                                orb_j1 = orb_is[ind_j]
                                orb_j2 = orb_js[ind_j]

                                J += ( scratch[ind_i,ind_j]
                                    * frame[i,orb_i1,k1] * frame[i,orb_i2,k2]
                                    * frame[i,orb_j1,j1] * frame[i,orb_j2,j2] )

                                if (orb_i1 != orb_i2):
                                    J += ( scratch[ind_i,ind_j]
                                        * frame[i,orb_i2,k1] * frame[i,orb_i1,k2]
                                        * frame[i,orb_j1,j1] * frame[i,orb_j2,j2] )

                                    if (orb_j1 != orb_j2):
                                        J += ( scratch[ind_i,ind_j]
                                            * frame[i,orb_i2,k1] * frame[i,orb_i1,k2]
                                            * frame[i,orb_j2,j1] * frame[i,orb_j1,j2] )

                                if (orb_j1 != orb_j2):
                                    J += ( scratch[ind_i,ind_j]
                                        * frame[i,orb_i1,k1] * frame[i,orb_i2,k2]
                                        * frame[i,orb_j2,j1] * frame[i,orb_j1,j2] )

                        epsilon_i[j1,j2,k1,k2] = J
                        epsilon_i[j2,j1,k1,k2] = J
                        epsilon_i[j1,j2,k2,k1] = J
                        epsilon_i[j2,j1,k2,k1] = J
    return epsilon

### jitted get F2e: the Fock matrix 'F'
@jit(nopython=True,parallel=True,fastmath=True)
def inline_F2e(F, P, epsilon, g, h, orbatm, orbnum, norbs, idxorbs):
    def tri_idx(i,j,N):
        return int( (N-1-(i+1)/2)*i + j - 0.5)

    natms = norbs.shape[0]
    norbstot = P.shape[0]

    for i in prange(norbstot):
        for j in prange(i,norbstot):
            atmi = orbatm[i]
            atmj = orbatm[j]
            Fij = 0
            if (atmi==atmj):
                if (i==j):
                    for k in range(idxorbs[atmi],idxorbs[atmi]+norbs[atmi]):
                        Fij += P[k,k] * (g[k,i] - 0.5*h[k,i])
                else:
                    Fij = 0.5 * P[i,j] * (3.0 * h[i,j] - g[i,j])

                for atmk in range(atmi):
                    idx_ik = tri_idx(atmk,atmi,natms)
                    # all atmi,atmk pairs
                    for k in range(idxorbs[atmk],idxorbs[atmk]+norbs[atmk]):
                        Fij += P[k,k] * epsilon[idx_ik, orbnum[i], orbnum[j], orbnum[k], orbnum[k]]
                        for l in range(k+1,idxorbs[atmk]+norbs[atmk]):
                            Fij += 2 * P[k,l] * epsilon[idx_ik, orbnum[i], orbnum[j], orbnum[k], orbnum[l]]

                for atmk in range(atmi+1,natms):
                    idx_ik = tri_idx(atmi,atmk,natms)
                    # all atmk,atmi pairs
                    for k in range(idxorbs[atmk],idxorbs[atmk]+norbs[atmk]):
                        Fij += P[k,k] * epsilon[idx_ik, orbnum[k], orbnum[k], orbnum[i], orbnum[j]]
                        for l in range(k+1,idxorbs[atmk]+norbs[atmk]):
                            Fij += 2 * P[k,l] * epsilon[idx_ik, orbnum[k], orbnum[l], orbnum[i], orbnum[j]]
            else: #atmi != atmj
                if (orbatm[j] < orbatm[i]):
                    idx_ij = tri_idx(orbatm[j],orbatm[i],natms)
                    for k in range(idxorbs[atmi],idxorbs[atmi]+norbs[atmi]):
                        for l in range(idxorbs[atmj],idxorbs[atmj]+norbs[atmj]):
                            Fij -= 0.5 * P[k,l] * epsilon[idx_ij, orbnum[i], orbnum[k], orbnum[j], orbnum[l]]
                else:
                    idx_ij = tri_idx(orbatm[i],orbatm[j],natms)
                    for k in range(idxorbs[atmj],idxorbs[atmj]+norbs[atmj]):
                        for l in range(idxorbs[atmi],idxorbs[atmi]+norbs[atmi]):
                            Fij -= 0.5 * P[k,l] * epsilon[idx_ij, orbnum[j], orbnum[k], orbnum[i], orbnum[l]]

            F[i,j] = F[j,i] = Fij
        #Fs[n] = Fij
        #return Fij
    #return Fs

### jitted E_core: core-core repulsion energy
@jit(nopython=True,fastmath=True)
def inline_Ecore(atnos, Qs, epsilon_s, alphas, Ks, Ls, Ms, Rs, idx_is, idx_js):
    Ec = 0.0

    for idx in range(len(Rs)):
        i = idx_is[idx]
        j = idx_js[idx]

        Rij = Rs[idx]
        epsij = epsilon_s[idx]

        Emndo = 0
        if (atnos[i] == 1 and (atnos[j] == 7 or atnos[j] == 8)):
            Emndo = epsij * (
                1 + numpy.exp(-alphas[i]*Rij) + Rij*numpy.exp(-alphas[j]*Rij) )
        elif (atnos[j] == 1 and (atnos[i] == 7 or atnos[i] == 8)):
            Emndo = epsij * (
                1 + Rij*numpy.exp(-alphas[i]*Rij) + numpy.exp(-alphas[j]*Rij) )
        else:
            Emndo = epsij * (
                1 + numpy.exp(-alphas[i]*Rij) + numpy.exp(-alphas[j]*Rij) )

        Eam1 =  1.0 / Rij * ((
            Ks[i,0]*numpy.exp(-Ls[i,0]*numpy.square((Rij-Ms[i,0]))) +
            Ks[i,1]*numpy.exp(-Ls[i,1]*numpy.square((Rij-Ms[i,1]))) +
            Ks[i,2]*numpy.exp(-Ls[i,2]*numpy.square((Rij-Ms[i,2]))) +
            Ks[i,3]*numpy.exp(-Ls[i,3]*numpy.square((Rij-Ms[i,3])))
        ) + (
            Ks[j,0]*numpy.exp(-Ls[j,0]*numpy.square((Rij-Ms[j,0]))) +
            Ks[j,1]*numpy.exp(-Ls[j,1]*numpy.square((Rij-Ms[j,1]))) +
            Ks[j,2]*numpy.exp(-Ls[j,2]*numpy.square((Rij-Ms[j,2]))) +
            Ks[j,3]*numpy.exp(-Ls[j,3]*numpy.square((Rij-Ms[j,3])))
        ))

        Ec += Qs[i] * Qs[j] * (Emndo + Eam1)
    return Ec

@jit(nopython=True,parallel=True,fastmath=True)
def inline_sumP(P, C):
    norbs,nmodes = C.shape
    for i in prange(norbs):
        for j in prange(norbs):
            P[i,j] = 0
            for k in range(nmodes):
                P[i,j] += 2.0 * C[i,k] * C[j,k]

# simple molecular representation
class AM1Molecule:
    def __init__(self, filestream='', streamtype="pdb"):
        self.atoms = []

        self.orbitals = []
        self.active_orbitals = []

        self.name = 'LG1'
        self.chnid = 'A'
        if (filestream!=''):
            if( streamtype == "pdb" ):
                self.read_pdb(filestream)
            elif ( streamtype == 'mol2' ):
                  self.read_mol2(filestream)

    def read_pdb(self, filestream):
        for line in filestream:
            if line[:6] == "HETATM":
                self.name = line[17:20]
                self.chnid = line[21:22]
                elt = self.elt_from_name(line[12:16])
                atno = AM1db['atno'][elt]
                zetaS = AM1db['zeta.s'][elt]
                zetaP = AM1db['zeta.p'][elt]
                fields = (
                    line[12:16], elt,
                    (line[30:38], line[38:46], line[46:54]),
                    line[60:66],
                    atno, zetaS, zetaP)
                self.atoms.append( fields )

        self.atoms = numpy.array( self.atoms, dtype=atom_dtype )
        self.calculate_orbitals()

    def read_mol2(self, filestream):
        read_cont = False
        for line in filestream:
            if line.startswith('@<TRIPOS>ATOM'):
                read_cont = True
                continue
            if line.startswith('@<TRIPOS>BOND'):
                break
            if not read_cont: continue

            self.name = line[17:20] #?
            elt = self.elt_from_name(line[7:11])
            atno = AM1db['atno'][elt]
            zetaS = AM1db['zeta.s'][elt]
            zetaP = AM1db['zeta.p'][elt]
            fields = (
                line[7:11].strip(), elt,
                (line[16:26], line[26:36], line[36:46]),
                "1.0", #?
                atno, zetaS, zetaP)
            self.atoms.append( fields )

        print("read mol2!", len(self.atoms))
        self.atoms = numpy.array( self.atoms, dtype=atom_dtype )
        self.calculate_orbitals()
                  
    def norbitals(self):
        nrb = [ AM1db['norb'][i] for i in self.atoms['elt'] ]
        return numpy.sum(nrb)

    def nelectrons(self):
        nelec = [ AM1db['Z'][i] for i in self.atoms['elt'] ]
        return numpy.sum(nelec)

    def write_pdb(self, output):
        for atom in self.atoms:
            output.write("%-6s%5d %4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00%6.3f" % (
                'HETATM',1,atom["name"],self.name,'A',1, atom["X"][0],atom["X"][1],atom["X"][2],atom["B"]))
            #file=output )

    def elt_from_name(self, name):
        # special case: 'CA' in col 2
        if (name[1:2] == 'CA'):
            return
        name=name.strip()
        if (name in AM1db['Z']):
            return name
        elif (name[:-1] in AM1db['Z']):
            return name[:-1]
        elif (name[:-2] in AM1db['Z']):
            return name[:-2]
        elif (name[:-3] in AM1db['Z']):
            return name[:-3]
        elif (name[1:] in AM1db['Z']):
            return name[1:]
        elif (name[1:-1] in AM1db['Z']):
            return name[1:-1]
        elif (name[1:-2] in AM1db['Z']):
            return name[1:-2]
        else:
            print ("ERROR, unrecognized! ",name, name[1:-2])
            return ''

    # lookup orbital properties and populate orbital arrays
    def calculate_orbitals(self):
        maxorb = AM1db['maxorb']
        orbital = []
        natoms = self.atoms.shape[0]

        self.orbitals = numpy.zeros( (natoms,maxorb), dtype=orbital_dtype );
        self.active_orbitals = numpy.zeros( (natoms,maxorb), dtype=bool );
        for i in range(natoms):
            elt = self.atoms[i]['elt']
            norb = AM1db['norb'][elt]
            for j in range(norb):
                U = 0.0
                alpha = AM1db['alpha'][elt]
                beta = 0.0
                if (j == 0):
                    beta = AM1db['beta.s'][elt]
                    U = AM1db['U.s'][elt]
                else:
                    beta = AM1db['beta.p'][elt]
                    U = AM1db['U.p'][elt]
                self.orbitals[i,j] = ( i, j, U, alpha, beta)
                self.active_orbitals[i,j] = True

class AM1:
    def __init__(self, molstream, streamtype="pdb"):
        self.mol = AM1Molecule(molstream, streamtype=streamtype)

        self.num_deriv_step = 0.0001
        self.history_size = 3

        # initialization
        # atoms
        self.natms = len(self.mol.atoms)

        # orbitals
        self.orblist = self.mol.orbitals[self.mol.active_orbitals]
        self.norbs = len(self.orblist)
        self.orbatm,self.orbnum = numpy.nonzero(self.mol.active_orbitals)
        # (cache indices)
        self.norbs_i = numpy.sum(self.mol.active_orbitals,axis=1)
        self.idxorbs_i = numpy.cumsum(self.norbs_i) - self.norbs_i

        self.maxorbitals = self.mol.orbitals.shape[1] # max orbitals/atom
        orb_js, orb_is = numpy.meshgrid( numpy.arange(self.maxorbitals), numpy.arange(self.maxorbitals) )
        self.orb_ijs = numpy.nonzero( orb_is<=orb_js )

        # atom-pairs
        self.j_atmmask, self.i_atmmask = numpy.nonzero(
            numpy.arange(self.natms)[:,None] < numpy.arange(self.natms)[None,:])
        self.natmpairs = len(self.i_atmmask)

        # atm-pair distances
        self.Rvecs = numpy.zeros((self.natms,self.natms,3))
        self.Rs = numpy.zeros((self.natms,self.natms))
        self.frame = numpy.zeros((self.natmpairs,4,4))

        # intermediate storage space
        self.epsilon = numpy.zeros((self.natmpairs,4,4,4,4))
        self.S = numpy.eye(self.norbs)
        self.F = numpy.zeros((self.norbs,self.norbs))
        self.P = numpy.zeros((self.norbs,self.norbs))
        self.H = numpy.zeros((self.norbs,self.norbs))
        self.Ec = 0.0

        # now lookup and store parameters that only depend on atom types
        self.lookup_gh()
        self.lookup_multipoles()

        self.K_is = numpy.zeros((self.natms,4))
        self.L_is = numpy.zeros((self.natms,4))
        self.M_is = numpy.zeros((self.natms,4))
        self.Q_is = numpy.zeros(self.natms)
        for i in range(self.natms):
            self.K_is[i,:] = AM1db['K'][self.mol.atoms[i]['elt']]
            self.L_is[i,:] = AM1db['L'][self.mol.atoms[i]['elt']]
            self.M_is[i,:] = AM1db['M'][self.mol.atoms[i]['elt']]
            self.Q_is[i] = AM1db['Z'][self.mol.atoms[i]['elt']]

        self.Us = numpy.zeros( (self.norbs) )
        for i in range(self.norbs):
            if (self.orblist[i]['type'] == 0):
                self.Us[i] = AM1db['U.s'][self.mol.atoms[self.orblist[i]['atomno']]['elt']]
            else:
                self.Us[i] = AM1db['U.p'][self.mol.atoms[self.orblist[i]['atomno']]['elt']]

        self.Zs = numpy.zeros( self.natms )
        for i in range(self.natms):
            self.Zs[i] = AM1db['Z'][self.mol.atoms[i]['elt']]

        self.Eatmheat=0
        for i in range(self.natms):
            self.Eatmheat += (
                AM1db['Eheat'][self.mol.atoms[i]['elt']] -
                EVtoKcal * AM1db['Eisol'][self.mol.atoms[i]['elt']] )

        self.E_hf = 0
        self.E_scf = 0

    def lookup_gh(self):
        self.g = numpy.zeros((self.norbs,self.norbs))
        self.h = numpy.zeros((self.norbs,self.norbs))

        selfmask = (self.orblist['atomno'] == self.orblist['atomno'][:,None])
        elts = self.mol.atoms[self.orblist['atomno']]['elt']
        type1 = self.orblist['type']
        type2 = self.orblist['type'][:,None]

        # ppp = same p-p
        pppmask = (type1>0) & (type2>0) & (type1!=type2) & selfmask
        # pp = diff p-p
        ppmask = (type1>0) & (type2>0) & selfmask & numpy.bitwise_not(pppmask)
        # ss
        ssmask = (type1==0) & (type2==0) & selfmask
        # sp
        spmask = ((type1==0) | (type2==0)) & selfmask & numpy.bitwise_not(ssmask)

        # lookup params
        for uniqelt in numpy.unique(elts):
            eltmask = (elts==uniqelt) & selfmask
            self.g[eltmask & pppmask] = AM1db['g'][uniqelt]['ppp']
            self.h[eltmask & pppmask] = 0.5 * (
                AM1db['g'][uniqelt]['pp'] - AM1db['g'][uniqelt]['ppp'] )
            self.g[eltmask & ppmask] = AM1db['g'][uniqelt]['pp']
            self.h[eltmask & ppmask] = AM1db['g'][uniqelt]['pp']
            self.g[eltmask & spmask] = AM1db['g'][uniqelt]['sp']
            self.h[eltmask & spmask] = AM1db['h'][uniqelt]['sp']
            self.g[eltmask & ssmask] = AM1db['g'][uniqelt]['ss']
            self.h[eltmask & ssmask] = AM1db['g'][uniqelt]['ss']

    def lookup_multipoles(self):
        # convert between orbital indices and (i,j)s
        def tri_idx(i,j,N=4):
            return int( (N-(i+1)/2)*i + j + 0.5)

        self.P1s = numpy.zeros( (self.natmpairs, 10, 4, 5) )
        self.P2s = numpy.zeros( (self.natmpairs, 10, 4, 5) )

        for i in range(self.natmpairs):
            elt1 = self.mol.atoms[ self.j_atmmask[i] ]['elt']
            mp1 = AM1db['multipoles'][elt1]
            for j in range(len(mp1)):
                (p1i,p1j) = mp1[j][0]
                npoles = len(mp1[j][1])
                idx = tri_idx(p1i,p1j)
                self.P1s[i,idx,0:npoles,:] = numpy.array(mp1[j][1])

            elt2 = self.mol.atoms[ self.i_atmmask[i] ]['elt']
            mp2 = AM1db['multipoles'][elt2]
            for j in range(len(mp2)):
                (p2i,p2j) = mp2[j][0]
                npoles = len(mp2[j][1])
                idx = tri_idx(p2i,p2j)
                self.P2s[i,idx,0:npoles,:] = numpy.array(mp2[j][1])

    # get R/Rvecs: interatomic vectors/dists
    def get_Rs(self, atomno=-1):
        if (atomno == -1):
            iatm,jatm = numpy.triu_indices(self.natms,1)
        else:
            iatm = numpy.full( self.natms-1, atomno )
            jatm = numpy.concatenate( [numpy.arange(atomno), numpy.arange(atomno+1,self.natms)] )

        inline_Rs( iatm, jatm, self.mol.atoms['X'], self.Rvecs, self.Rs, self.frame )

    # calculates S: the orbital overlap matrix
    def get_S(self, atomno=-1):
        ovr.inline_S(
            atomno, self.S, self.mol.atoms['atno'],
            self.mol.atoms['zetaS'], self.mol.atoms['zetaP'],
            self.Rvecs, self.Rs, self.norbs_i, self.idxorbs_i )

    # Calculate epsilon, the two-center, two-electron integrals
    def get_epsilon(self, atomno=-1):
        if atomno != -1:
            atompairs, = numpy.where( (self.j_atmmask==atomno) | (self.i_atmmask==atomno) )
        else:
            atompairs = numpy.arange( self.natmpairs )

        self.epsilon[ atompairs ] = inline_epsilon(
            atompairs, self.P1s, self.P2s, self.orb_ijs[0], self.orb_ijs[1],
            self.j_atmmask, self.i_atmmask, self.Rs, self.Rvecs,
            self.frame
        )

    # Get initial Fock matrix using extended Huckel
    def get_F0(self):
        self.F = ((0.875)*(self.Us[:,None]+self.Us[None,:])) * self.S
        numpy.fill_diagonal(self.F, self.Us)
        return self.F

    # update F (the Fock matrix) from P (the density matrix)
    def update_F(self, P):
        inline_F2e(self.F, P, self.epsilon, self.g, self.h, self.orbatm,
                self.orbnum, self.norbs_i, self.idxorbs_i)
        return self.F

    # Calculate H: the one-electron matrix
    def get_H(self, atomno=-1):
        bs = self.orblist['beta']
        self.H = 0.5 * (bs[:,None] + bs[None,:]) * self.S

        # fill in single-center elts
        inline_V(
            self.H, self.Zs, self.epsilon,
            self.orbatm, self.orbnum, self.norbs_i, self.idxorbs_i)

        # append diag
        self.H += numpy.diag(self.Us)

        return

    # Calculate P: the density matrix
    def get_P(self, F):
        Fp = F + self.H
        _, C = numpy.linalg.eigh(Fp)
        npairs = int(self.mol.nelectrons() / 2) # verify even?
        inline_sumP( self.P, C[:,:npairs] )
        return self.P

    # Calculate the core-core repulsion
    def get_Ecore(self):
        myRs = self.Rs[self.i_atmmask,self.j_atmmask]*BohrRadius
        self.Ec = inline_Ecore(
            self.mol.atoms['atno'],
            self.Q_is, self.epsilon[:,0,0,0,0], self.mol.orbitals[:,0]['alpha'],
            self.K_is, self.L_is, self.M_is, myRs, self.i_atmmask, self.j_atmmask
        )
        return self.Ec

    def update_Energy(self, F, P):
        Eel = numpy.sum( 0.5*P * (self.H + F) )
        E = (Eel + self.Ec)
        self.E_scf = E*EVtoKcal
        self.E_hf = E*EVtoKcal + self.Eatmheat
        return E

    # Calculates the energy by doing a full SCF
    def get_Energy(self):
        E1 = 0.0

        self.get_Rs()
        self.get_S()
        self.get_epsilon()
        self.get_H()
        F0 = self.get_F0()
        P1 = self.get_P(F0)
        self.get_Ecore()
        G = self.update_F(P1)
        E2 = self.update_Energy(G + self.H, P1)

        (damp, E_down, maxiter, tol) = (0., 0, 400, 1e-6)
        niter = 0
        while (True):
            niter += 1
            P0 = P1.copy()
            P1 = self.get_P(G)
            P1 = (1 - damp) * P1 + damp * P0
            G = self.update_F(P1)
            E1 = E2
            E2 = self.update_Energy(G + self.H, P1)

            if E2 > E1:
                damp = 0.5 * (1 + damp)
                E_down = 0
            else:
                E_down += 1
                if E_down > 3:
                    damp *= 0.5

            if (niter >= maxiter or numpy.max(numpy.abs(P1 - P0)) <= tol):
                break

        self.G = G
        self.P = P1
        if niter >= maxiter:
            print('WARNING: SCF reached iter %i :'% maxiter, E2 , E1 )
        return E2

    # Resets energy during minimization
    def reset_Energy_ext(self, P, atomno=-1):
        self.get_Rs(atomno)
        self.get_S(atomno)
        self.get_epsilon(atomno)

    # Calculates energy during minimization
    def get_Energy_ext(self, atomno=-1):
        self.get_Rs(atomno)
        self.get_S(atomno)
        self.get_epsilon(atomno)
        self.get_H()
        self.get_Ecore()
        G = self.update_F(self.P)
        E = self.update_Energy(G + self.H, self.P)
        return E

    # get Mullikan charges
    def get_charges(self):
        natms = len(self.mol.atoms)
        orblist = self.mol.orbitals[self.mol.active_orbitals]
        self.Zs = numpy.empty(natms)
        for i in range(natms):
            self.Zs[i] = (
                AM1db['Z'][ self.mol.atoms[i]['elt'] ]
                - numpy.sum( numpy.diag(self.P)[orblist['atomno']==i] )
            )
        return self.Zs


    # gradient (numeric)
    def gradient(self):
        natms = len(self.mol.atoms)

        E0 = self.get_Energy()
        Rref = numpy.copy( self.Rvecs )
        grad = numpy.zeros( (natms,3) )

        for i in range(natms):
            Xi0 = numpy.copy( self.mol.atoms[i]['X'] )
            for j in range(3):
                self.mol.atoms[i]['X'][j] = Xi0[j] + self.num_deriv_step
                Eplus = self.get_Energy_ext(i)
                self.mol.atoms[i]['X'][j] = Xi0[j] - self.num_deriv_step
                Eminus = self.get_Energy_ext(i)
                self.mol.atoms[i]['X'][j] = Xi0[j]
                self.reset_Energy_ext(self.P, i) # reset coords

                grad[i,j] = (Eplus - Eminus) / (2*self.num_deriv_step)

        return grad

    def armijo_linesearch(
            self, func, E0, dir_deriv, alpha0=1,
            factor=0.5, sigma_decrease=1e-4, sigma_increase=0.5, minstep=1e-18
    ):
        E_a0 = func(alpha0)

        if E_a0 <= E0 + alpha0 * sigma_increase * dir_deriv:
            # attempt to increase stepsize
            alpha1 = alpha0 / factor
            E_a1 = func(alpha1)
            if (E_a1 < E_a0):
                return E_a1, alpha1

            # step back
            return E_a0, alpha0

        # next, we check if we need to decrease the stepsize
        alpha1 = alpha0
        E_a1 = E_a0
        while (E_a1 > E0 + alpha1 * sigma_decrease * dir_deriv):
            if (alpha1 < minstep):
                if (E_a1 >= E0):
                    finite_diff = (E_a1 - E0) / alpha1
                    print(
                        "Inaccurate G! Step=", alpha1, " Deriv=", dir_deriv,
                        " Finite=", finite_diff
                    )
                    return E0, 0.0
                return alpha1, E_a1
            alpha1 *= factor
            E_a1 = func(alpha1)

        return E_a1, alpha1

    # l-bfgs
    def runLBFGS( self, maxiter=1000, tol_g=1e-4, tol_e=1e-8 ):
        shift = self.mol.atoms['X'][0,:]

        def getDofs():
            return self.mol.atoms['X'].reshape(-1)
        def setDofs(X):
            self.mol.atoms['X'] = X.reshape(-1,3)
        def getdFunc():
            return self.gradient().reshape(-1)

        niter = 0
        loss = self.get_Energy()
        norm_g = 2*tol_g
        laststep = 2*tol_e
        dofs = getDofs()

        g = numpy.zeros_like( dofs )

        while niter < maxiter:
            niter += 1

            g0 = g
            g = getdFunc()
            max_grad = g.max()

            # converge check 1: gradient
            if max_grad <= tol_g:
                #print ("Converged with max(g) = ",max_grad)
                break

            if niter == 1:
                d = -g
                old_dirs = []
                old_stps = []
                ro = [None] * self.history_size
                al = [None] * self.history_size
                t = 1.0 / numpy.sqrt( numpy.dot(d,d) )
            else:
                # lbfgs update (update memory)
                y = g - g0
                s = d*t
                ys = y.dot(s)
                yy = y.dot(y)
                if ys > 1e-10:
                    if len(old_dirs) == self.history_size:
                        old_dirs.pop(0)
                        old_stps.pop(0)
                    old_dirs.append(y)
                    old_stps.append(s)

                num_old = len(old_dirs)

                for i in range(num_old):
                    ro[i] = 1. / old_dirs[i].dot(old_stps[i])

                # iteration in L-BFGS loop collapsed to use just one buffer
                d = -g
                for i in range(num_old - 1, -1, -1):
                    al[i] = old_stps[i].dot(d) * ro[i]
                    d += -al[i]*old_dirs[i]

                d *= (ys/yy)

                for i in range(num_old):
                    be_i = old_dirs[i].dot(d) * ro[i]
                    d += (al[i] - be_i)*old_stps[i]

            prev_loss = loss
            gtd = g.dot(d)  # g * d

            def linefunc(alpha):
                setDofs(alpha*d + dofs)
                E = self.get_Energy()
                setDofs(dofs)
                return E

            loss, t = self.armijo_linesearch(linefunc, prev_loss, gtd, t)
            setDofs(t*d + dofs)
            dofs = getDofs()
            if (niter % 10 == 0):
                print ("iter ",niter," E=",loss*EVtoKcal,"  step=",t)

            # converge check 2: rel tol
            if 2 * abs(loss - prev_loss) <= tol_e * (
                    abs(loss) + abs(prev_loss) + 1e-10):
                break

        print ("final E=",loss*EVtoKcal,"  step=",t)
        return loss


## CH4
# Heat of formation   =         -8.79856886 kcal/mol  (       -0.38153458 eV)
# Total SCF energy    =      -4225.48160648 kcal/mol  (     -183.23063208 eV)
# Electronic energy   =      -8844.18839418 kcal/mol  (     -383.51278757 eV)
# Core-core repulsion =       4618.70678770 kcal/mol  (      200.28215549 eV)
def testCH4():
    molstring = (
        "HETATM    1  C1  CH4     1       0.000   0.000   0.000  1.00  0.00\n"
        "HETATM    2  H1  CH4     1       1.000   0.000   0.000  1.00  0.00\n"
        "HETATM    3  H2  CH4     1       1.000   1.000   0.000  1.00  0.00\n"
        "HETATM    4  H3  CH4     1      -1.000   0.000   1.000  1.00  0.00\n"
        "HETATM    5  H4  CH4     1      -1.000   0.000   0.000  1.00  0.00\n"
    );
    molstream = io.StringIO(molstring)
    Am1Mol = AM1(molstream)
    Am1Mol.runLBFGS()
    print ("Heat of formation = ", Am1Mol.E_hf, " kcal/mol")
    print ("Total SCF energy = ", Am1Mol.E_scf, " kcal/mol")
    print ("Electronic energy = ", Am1Mol.E_scf-Am1Mol.Ec*EVtoKcal, " kcal/mol")
    print ("Core-core repulsion = ", Am1Mol.Ec*EVtoKcal, " kcal/mol")


## NH3
# Heat of formation   =         -7.29251288 kcal/mol  (       -0.31622709 eV)
# Total SCF energy    =      -5732.76248335 kcal/mol  (     -248.59123556 eV)
# Electronic energy   =     -10024.63544328 kcal/mol  (     -434.70081277 eV)
# Core-core repulsion =       4291.87295993 kcal/mol  (      186.10957721 eV)
def testNH3():
    molstring = (
        "HETATM    1  N1  NH3     1       0.103   0.103   0.103  1.00  0.00\n"
        "HETATM    2  H1  NH3     1       1.065  -0.083  -0.085  1.00  0.00\n"
        "HETATM    3  H2  NH3     1      -0.085   1.065  -0.083  1.00  0.00\n"
        "HETATM    4  H3  NH3     1      -0.083  -0.086   1.065  1.00  0.00\n"
    );
    molstream = io.StringIO(molstring)
    Am1Mol = AM1(molstream)
    Am1Mol.runLBFGS()
    #Am1Mol.get_Energy()
    #Am1Mol.mol.write_pdb(None)
    print ("Heat of formation = ", Am1Mol.E_hf, " kcal/mol")
    print ("Total SCF energy = ", Am1Mol.E_scf, " kcal/mol")
    print ("Electronic energy = ", Am1Mol.E_scf-Am1Mol.Ec*EVtoKcal, " kcal/mol")
    print ("Core-core repulsion = ", Am1Mol.Ec*EVtoKcal, " kcal/mol")

## CNH
# Heat of formation   =         30.98097069 kcal/mol  (        1.34343570 eV)
# Total SCF energy    =      -8021.68201912 kcal/mol  (     -347.84623473 eV)
# Electronic energy   =     -14095.11241184 kcal/mol  (     -611.20993937 eV)
# Core-core repulsion =       6073.43039272 kcal/mol  (      263.36370464 eV)
def testCNH():
    molstring = (
        "HETATM    1  C1  HCN     1      0.0330  0.000   0.000  1.00  0.00\n"
        "HETATM    2  N2  HCN     1      1.1930  0.000   0.000  1.00  0.00\n"
        "HETATM    3  H1  HCN     1     -1.0360  0.000   0.000  1.00  0.00\n"
    );
    molstream = io.StringIO(molstring)
    Am1Mol = AM1(molstream)
    Am1Mol.runLBFGS()
    #Am1Mol.get_Energy()
    print ("Heat of formation = ", Am1Mol.E_hf, " kcal/mol")
    print ("Total SCF energy = ", Am1Mol.E_scf, " kcal/mol")
    print ("Electronic energy = ", Am1Mol.E_scf-Am1Mol.Ec*EVtoKcal, " kcal/mol")
    print ("Core-core repulsion = ", Am1Mol.Ec*EVtoKcal, " kcal/mol")

## CO2
# Heat of formation   =        -79.85232260 kcal/mol  (       -3.46265655 eV)
# Total SCF energy    =     -17735.13540947 kcal/mol  (     -769.05318111 eV)
# Electronic energy   =     -33118.43934618 kcal/mol  (    -1436.12329674 eV)
# Core-core repulsion =      15383.30393671 kcal/mol  (      667.07011564 eV)
def testCO2():
    molstring = (
        "HETATM    1  C1  CO2     1       0.000   0.000   0.000  1.00  0.00\n"
        "HETATM    2  O1  CO2     1       1.300   0.000   0.000  1.00  0.00\n"
        "HETATM    2  O1  CO2     1      -1.300   0.000   0.000  1.00  0.00\n"
    );
    molstream = io.StringIO(molstring)
    Am1Mol = AM1(molstream)
    Am1Mol.runLBFGS()
    print ("Heat of formation = ", Am1Mol.E_hf, " kcal/mol")
    print ("Total SCF energy = ", Am1Mol.E_scf, " kcal/mol")
    print ("Electronic energy = ", Am1Mol.E_scf-Am1Mol.Ec*EVtoKcal, " kcal/mol")
    print ("Core-core repulsion = ", Am1Mol.Ec*EVtoKcal, " kcal/mol")

## Alanine
# Heat of formation   =        -66.32493822 kcal/mol  (       -2.87606514 eV)
# Total SCF energy    =     -30620.78405862 kcal/mol  (    -1327.81683616 eV)
# Electronic energy   =    -107272.62299942 kcal/mol  (    -4651.68999607 eV)
# Core-core repulsion =      76651.83894080 kcal/mol  (     3323.87315992 eV)
def testALA():
    molstring = (
        "HETATM    1  N   ALA X   1       0.227  -1.259   0.452  1.00  0.50\n"
        "HETATM    2  CA  ALA X   1       0.103  -0.030  -0.392  1.00  0.50\n"
        "HETATM    3  C   ALA X   1       1.270   0.922  -0.094  1.00  0.50\n"
        "HETATM    4  O   ALA X   1       2.008   1.323  -0.994  1.00  0.50\n"
        "HETATM    5  CB  ALA X   1      -1.244   0.625  -0.159  1.00  0.50\n"
        "HETATM    6  OXT ALA X   1       1.498   1.305   1.054  1.00  0.50\n"
        "HETATM    7  H   ALA X   1       0.069  -1.019   1.421  1.00  0.50\n"
        "HETATM    8  HA  ALA X   1       0.160  -0.299  -1.339  1.00  0.50\n"
        "HETATM    9  HB1 ALA X   1      -1.150   1.385   0.442  1.00  0.50\n"
        "HETATM   10  HB2 ALA X   1      -1.605   0.932  -1.008  1.00  0.50\n"
        "HETATM   11  HB3 ALA X   1      -1.857  -0.018   0.234  1.00  0.50\n"
        "HETATM   12  H2  ALA X   1       1.104  -1.640   0.356  1.00  0.50\n"
        "HETATM   13  H3  ALA X   1      -0.424  -1.909   0.174  1.00  0.50\n"
    );
    molstream = io.StringIO(molstring)
    Am1Mol = AM1(molstream)
    Am1Mol.runLBFGS()
    Z = Am1Mol.get_charges()
    Am1Mol.mol.atoms['B'] = Z
    Am1Mol.mol.write_pdb(open("ala_out.pdb", "w"))
    print ("Heat of formation = ", Am1Mol.E_hf, " kcal/mol")
    print ("Total SCF energy = ", Am1Mol.E_scf, " kcal/mol")
    print ("Electronic energy = ", Am1Mol.E_scf-Am1Mol.Ec*EVtoKcal, " kcal/mol")
    print ("Core-core repulsion = ", Am1Mol.Ec*EVtoKcal, " kcal/mol")

## Tryptophan
# Heat of formation   =        -12.10828652 kcal/mol  (       -0.52505470 eV)
# Total SCF energy    =     -60578.05158695 kcal/mol  (    -2626.86143649 eV)
# Electronic energy   =    -339934.29238393 kcal/mol  (   -14740.65705667 eV)
# Core-core repulsion =     279356.24079698 kcal/mol  (    12113.79562018 eV)
def testTRP():
    molstring = (
        "HETATM    1  N   TRP X   1       0.224  -1.242   0.445  1.00  0.00\n"
        "HETATM    2  CA  TRP X   1       0.129  -0.052  -0.396  1.00  0.00\n"
        "HETATM    3  C   TRP X   1       1.247   0.917  -0.092  1.00  0.00\n"
        "HETATM    4  O   TRP X   1       2.008   1.323  -0.994  1.00  0.50\n"
        "HETATM    5  CB  TRP X   1      -1.243   0.635  -0.151  1.00  0.00\n"
        "HETATM    6  CG  TRP X   1      -2.452  -0.135  -0.691  1.00  0.00\n"
        "HETATM    7  CD1 TRP X   1      -3.031   0.032  -1.966  1.00  0.00\n"
        "HETATM    8  CD2 TRP X   1      -3.248  -1.048  -0.031  1.00  0.00\n"
        "HETATM    9  CE2 TRP X   1      -4.294  -1.427  -0.909  1.00  0.00\n"
        "HETATM   10  CE3 TRP X   1      -3.181  -1.579   1.283  1.00  0.00\n"
        "HETATM   11  NE1 TRP X   1      -4.179  -0.768  -2.124  1.00  0.00\n"
        "HETATM   12  CZ2 TRP X   1      -5.274  -2.350  -0.483  1.00  0.00\n"
        "HETATM   13  CZ3 TRP X   1      -4.157  -2.494   1.678  1.00  0.00\n"
        "HETATM   14  CH2 TRP X   1      -5.188  -2.876   0.808  1.00  0.00\n"
        "HETATM   15  OXT TRP X   1       1.498   1.305   1.054  1.00  0.50\n"
        "HETATM   16  H   TRP X   1       0.069  -1.019   1.421  1.00  0.50\n"
        "HETATM   17  HA  TRP X   1       0.224  -0.350  -1.457  1.00  0.00\n"
        "HETATM   18 2HB  TRP X   1      -1.256   1.641  -0.613  1.00  0.00\n"
        "HETATM   19 3HB  TRP X   1      -1.385   0.829   0.931  1.00  0.00\n"
        "HETATM   20 1HD  TRP X   1      -2.662   0.727  -2.709  1.00  0.00\n"
        "HETATM   21 1HE  TRP X   1      -4.815  -0.816  -2.928  1.00  0.00\n"
        "HETATM   22 3HE  TRP X   1      -2.394  -1.284   1.963  1.00  0.00\n"
        "HETATM   23 2HZ  TRP X   1      -6.074  -2.641  -1.148  1.00  0.00\n"
        "HETATM   24 3HZ  TRP X   1      -4.116  -2.914   2.673  1.00  0.00\n"
        "HETATM   25 2HH  TRP X   1      -5.926  -3.589   1.141  1.00  0.00\n"
        "HETATM   26  H2  TRP X   1       1.104  -1.640   0.356  1.00  0.50\n"
        "HETATM   27  H3  TRP X   1      -0.424  -1.909   0.174  1.00  0.50\n"
    );
    molstream = io.StringIO(molstring)
    Am1Mol = AM1(molstream)
    Am1Mol.runLBFGS()
    Z = Am1Mol.get_charges()
    Am1Mol.mol.atoms['B'] = Z
    Am1Mol.mol.write_pdb(open("trp_out.pdb", "w"))
    print ("Heat of formation = ", Am1Mol.E_hf, " kcal/mol")
    print ("Total SCF energy = ", Am1Mol.E_scf, " kcal/mol")
    print ("Electronic energy = ", Am1Mol.E_scf-Am1Mol.Ec*EVtoKcal, " kcal/mol")
    print ("Core-core repulsion = ", Am1Mol.Ec*EVtoKcal, " kcal/mol")

## Big case (plus needs significant movement)
# Heat of formation   =       -241.89202266 kcal/mol  (      -10.48922521 eV)
# Total SCF energy    =    -155471.24297255 kcal/mol  (    -6741.73899538 eV)
# Electronic energy   =   -1591627.52331202 kcal/mol  (   -69018.14853268 eV)
# Core-core repulsion =    1436156.28033947 kcal/mol  (    62276.40953729 eV)
def testBIG():
    molstring = (
        "HETATM    1  C1  UNK     1      -9.244   0.248   2.613  1.00 20.00\n"
        "HETATM    2  C2  UNK     1      -7.765   0.630   2.523  1.00 20.00\n"
        "HETATM    3  C3  UNK     1      -7.109  -0.139   1.374  1.00 20.00\n"
        "HETATM    4  C4  UNK     1      -5.630   0.243   1.284  1.00 20.00\n"
        "HETATM    5  N5  UNK     1      -5.002  -0.493   0.184  1.00 20.00\n"
        "HETATM    6  C6  UNK     1      -3.694  -0.316  -0.086  1.00 20.00\n"
        "HETATM    7  O7  UNK     1      -3.036   0.454   0.583  1.00 20.00\n"
        "HETATM    8  C8  UNK     1      -3.048  -1.073  -1.217  1.00 20.00\n"
        "HETATM    9  H9  UNK     1      -3.546  -0.821  -2.154  1.00 20.00\n"
        "HETATM   10  C10 UNK     1      -3.171  -2.576  -0.961  1.00 20.00\n"
        "HETATM   11  C11 UNK     1      -1.569  -0.691  -1.308  1.00 20.00\n"
        "HETATM   12  C12 UNK     1      -1.445   0.774  -1.730  1.00 20.00\n"
        "HETATM   13  H13 UNK     1      -2.033   1.398  -1.056  1.00 20.00\n"
        "HETATM   14  O14 UNK     1      -1.931   0.928  -3.065  1.00 20.00\n"
        "HETATM   15  C15 UNK     1       0.024   1.201  -1.669  1.00 20.00\n"
        "HETATM   16  H16 UNK     1       0.397   1.092  -0.655  1.00 20.00\n"
        "HETATM   17  C17 UNK     1       0.146   2.662  -2.112  1.00 20.00\n"
        "HETATM   18  C18 UNK     1      -0.142   3.600  -0.940  1.00 20.00\n"
        "HETATM   19  H19 UNK     1      -0.373   3.013  -0.053  1.00 20.00\n"
        "HETATM   20  C20 UNK     1      -1.336   4.494  -1.284  1.00 20.00\n"
        "HETATM   21  C21 UNK     1       1.077   4.481  -0.662  1.00 20.00\n"
        "HETATM   22  C22 UNK     1       1.532   4.274   0.785  1.00 20.00\n"
        "HETATM   23  C23 UNK     1       2.505   3.094   0.837  1.00 20.00\n"
        "HETATM   24  C24 UNK     1       3.770   3.484   1.600  1.00 20.00\n"
        "HETATM   25  C25 UNK     1       4.332   2.258   2.326  1.00 20.00\n"
        "HETATM   26  N26 UNK     1       4.197   1.071   1.452  1.00 20.00\n"
        "HETATM   27  C27 UNK     1       4.904   0.994   0.167  1.00 20.00\n"
        "HETATM   28  C28 UNK     1       6.278   0.355   0.378  1.00 20.00\n"
        "HETATM   29  C29 UNK     1       3.441   0.048   1.930  1.00 20.00\n"
        "HETATM   30  O30 UNK     1       3.253  -0.064   3.135  1.00 20.00\n"
        "HETATM   31  C31 UNK     1       2.787  -0.904   1.010  1.00 20.00\n"
        "HETATM   32  C32 UNK     1       2.735  -2.261   1.354  1.00 20.00\n"
        "HETATM   33  C33 UNK     1       2.130  -3.171   0.503  1.00 20.00\n"
        "HETATM   34  C34 UNK     1       1.567  -2.752  -0.691  1.00 20.00\n"
        "HETATM   35  C35 UNK     1       1.613  -1.399  -1.044  1.00 20.00\n"
        "HETATM   36  C36 UNK     1       2.231  -0.486  -0.193  1.00 20.00\n"
        "HETATM   37  C37 UNK     1       0.995  -0.947  -2.306  1.00 20.00\n"
        "HETATM   38  O38 UNK     1       0.650  -1.785  -3.131  1.00 20.00\n"
        "HETATM   39  N39 UNK     1       0.802   0.354  -2.581  1.00 20.00\n"
        "HETATM   40  O40 UNK     1       2.089  -4.488   0.839  1.00 20.00\n"
        "HETATM   41  C41 UNK     1       1.454  -5.374  -0.084  1.00 20.00\n"
        "HETATM   42  C42 UNK     1       1.505  -6.803   0.461  1.00 20.00\n"
        "HETATM   43  H43 UNK     1      -9.742   0.500   1.677  1.00 20.00\n"
        "HETATM   44  H44 UNK     1      -9.332  -0.823   2.796  1.00 20.00\n"
        "HETATM   45  H45 UNK     1      -9.711   0.796   3.431  1.00 20.00\n"
        "HETATM   46  H46 UNK     1      -7.267   0.379   3.459  1.00 20.00\n"
        "HETATM   47  H47 UNK     1      -7.677   1.701   2.340  1.00 20.00\n"
        "HETATM   48  H48 UNK     1      -7.607   0.113   0.437  1.00 20.00\n"
        "HETATM   49  H49 UNK     1      -7.197  -1.210   1.556  1.00 20.00\n"
        "HETATM   50  H50 UNK     1      -5.132  -0.008   2.220  1.00 20.00\n"
        "HETATM   51  H51 UNK     1      -5.542   1.314   1.101  1.00 20.00\n"
        "HETATM   52  H52 UNK     1      -5.528  -1.108  -0.351  1.00 20.00\n"
        "HETATM   53  H53 UNK     1      -2.637  -3.122  -1.738  1.00 20.00\n"
        "HETATM   54  H54 UNK     1      -2.742  -2.816   0.012  1.00 20.00\n"
        "HETATM   55  H55 UNK     1      -4.223  -2.861  -0.974  1.00 20.00\n"
        "HETATM   56  H56 UNK     1      -1.074  -1.324  -2.044  1.00 20.00\n"
        "HETATM   57  H57 UNK     1      -1.097  -0.830  -0.334  1.00 20.00\n"
        "HETATM   58  H58 UNK     1      -1.453   0.401  -3.720  1.00 20.00\n"
        "HETATM   59  H59 UNK     1       1.154   2.838  -2.487  1.00 20.00\n"
        "HETATM   60  H60 UNK     1      -0.569   2.848  -2.913  1.00 20.00\n"
        "HETATM   61  H61 UNK     1      -2.211   3.874  -1.477  1.00 20.00\n"
        "HETATM   62  H62 UNK     1      -1.542   5.163  -0.448  1.00 20.00\n"
        "HETATM   63  H63 UNK     1      -1.105   5.083  -2.172  1.00 20.00\n"
        "HETATM   64  H64 UNK     1       1.888   4.202  -1.337  1.00 20.00\n"
        "HETATM   65  H65 UNK     1       0.818   5.528  -0.823  1.00 20.00\n"
        "HETATM   66  H66 UNK     1       2.021   5.178   1.144  1.00 20.00\n"
        "HETATM   67  H67 UNK     1       0.660   4.048   1.402  1.00 20.00\n"
        "HETATM   68  H68 UNK     1       2.021   2.247   1.320  1.00 20.00\n"
        "HETATM   69  H69 UNK     1       2.774   2.819  -0.186  1.00 20.00\n"
        "HETATM   70  H70 UNK     1       4.516   3.871   0.906  1.00 20.00\n"
        "HETATM   71  H71 UNK     1       3.523   4.256   2.336  1.00 20.00\n"
        "HETATM   72  H72 UNK     1       5.389   2.417   2.553  1.00 20.00\n"
        "HETATM   73  H73 UNK     1       3.780   2.093   3.253  1.00 20.00\n"
        "HETATM   74  H74 UNK     1       4.325   0.389  -0.530  1.00 20.00\n"
        "HETATM   75  H75 UNK     1       5.028   1.998  -0.240  1.00 20.00\n"
        "HETATM   76  H76 UNK     1       6.802   0.298  -0.576  1.00 20.00\n"
        "HETATM   77  H77 UNK     1       6.857   0.961   1.075  1.00 20.00\n"
        "HETATM   78  H78 UNK     1       6.153  -0.648   0.785  1.00 20.00\n"
        "HETATM   79  H79 UNK     1       3.171  -2.598   2.283  1.00 20.00\n"
        "HETATM   80  H80 UNK     1       1.092  -3.467  -1.344  1.00 20.00\n"
        "HETATM   81  H81 UNK     1       2.273   0.557  -0.477  1.00 20.00\n"
        "HETATM   82  H82 UNK     1       1.128   0.742  -3.404  1.00 20.00\n"
        "HETATM   83  H83 UNK     1       0.414  -5.074  -0.219  1.00 20.00\n"
        "HETATM   84  H84 UNK     1       1.971  -5.332  -1.043  1.00 20.00\n"
        "HETATM   85  H85 UNK     1       2.545  -7.103   0.596  1.00 20.00\n"
        "HETATM   86  H86 UNK     1       0.988  -6.845   1.420  1.00 20.00\n"
        "HETATM   87  H87 UNK     1       1.021  -7.479  -0.243  1.00 20.00\n"
    );
    molstream = io.StringIO(molstring)
    Am1Mol = AM1(molstream)
    Am1Mol.runLBFGS()
    Z = Am1Mol.get_charges()
    Am1Mol.mol.atoms['B'] = Z
    Am1Mol.mol.write_pdb(open("big_out.pdb", "w"))
    print ("Heat of formation = ", Am1Mol.E_hf, " kcal/mol")
    print ("Total SCF energy = ", Am1Mol.E_scf, " kcal/mol")
    print ("Electronic energy = ", Am1Mol.E_scf-Am1Mol.Ec*EVtoKcal, " kcal/mol")
    print ("Core-core repulsion = ", Am1Mol.Ec*EVtoKcal, " kcal/mol")

if __name__ == "__main__":
    #print("Testing CH4")
    #testCH4()
    #print("Testing NH3")
    #testNH3()
    #print("Testing CNH")
    #testCNH()
    #print("Testing CO2")
    #testCO2()
    #print("Testing Alanine")
    #testALA()
    print("Testing Tryptophan")
    testTRP()
    #print("Testing BIG")
    #testBIG()
