#!/usr/bin/env python3
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
'''
Defines a class assigning generic atom types into input Molecule class.
Uses hybridization, ring, bond, element type info from Molecule.

Author: Hahnbeom Park and Frank DiMaio 
'''

import sys

from Types import *
from BasicClasses import FunctionalGroupClass

class AtomTypeClassifier:
    def __init__(self):
        self.found_error = False
        return

    def apply_to_molecule(self,molecule):
        for iatm,atm in enumerate(molecule.atms):
            aclass = self.apply_atom(molecule,iatm) #integer
            molecule.atms[iatm].aclass = aclass

        # edit polar C type for special functional groups
        for iatm,atm in enumerate(molecule.atms):
            self.modify_polarC(molecule,iatm)

        #Correct amide bonds as border=4; this keep happens with AM1BCC...
        for i,bond in enumerate(molecule.bonds):
            aclass1 = ACLASS_ID[molecule.atms[bond.atm1].aclass]
            aclass2 = ACLASS_ID[molecule.atms[bond.atm2].aclass]
            is_amide = (aclass1 in ['Nad','Nad3'] and aclass2 == 'CDp') or \
                       (aclass2 in ['Nad','Nad3'] and aclass1 == 'CDp')
            if is_amide and bond.order == 1:
                #print('Warning: Amide bond b/w %s-%s was defined as single-bonded; correcting as amide-bond'%(aname1,aname2))
                molecule.bonds[i].order = 4

    def apply_atom(self,molecule,iatm): #input as MoleculeType
        atm = molecule.atms[iatm]
        # atype is string
        if atm.atype == 1: #C
            atype = self.classify_C(molecule,iatm)
        elif atm.atype == 2: #H
            atype = self.classify_H(molecule,iatm)
        elif atm.atype == 3: #O
            atype = self.classify_O(molecule,iatm)
            #if atype == 'OG2': print molecule.atms[iatm].name

        elif atm.atype == 4: #N
            atype = self.classify_N(molecule,iatm)

        elif atm.atype == 5: #P
            atype = self.classify_P(molecule,iatm)

        elif atm.atype == 6: #S
            atype = self.classify_S(molecule,iatm)

        elif atm.atype in [7,8,9,10]: #Halogen
            atype = self.classify_halogen(molecule,iatm)
        else:
            atype = self.classify_else(molecule,iatm)

        # return index
        if not atype:
            self.report_error(molecule,iatm)
            self.found_error = True
            return 0 #==Null
        else:
            return ACLASS_ID.index(atype)

    def report_error(self,molecule,iatm):
        print('ERROR AtomTyping: %10s %3d %4s'%(molecule.name,iatm,
                                                molecule.atms[iatm].name))
        
    # Defines rules for classification
    def classify_H(self,molecule,iatm):
        atm = molecule.atms[iatm]
        if len(atm.bonds) > 1:
            print('Warning, Hydrogen cannot have more than bond!')
            sys.exit()
            
        jatm,order = atm.bonds[0]
        jtype = molecule.atms[jatm].atype
        if jtype == 1: #C
            #backtrace
            #if molecule.atms[jatm].connected_to_polar:
            #    return 'HCW'
            if molecule.atms[jatm].hyb == 9:
                return 'HR'
            else:
                return 'HC'
        elif jtype == 3: #O
            return 'HO'
        elif jtype == 4: #N
            return 'HN'
        elif jtype == 6: #S
            return 'HS'
        else: 
            return 'HG'

    def classify_C(self,molecule,iatm):
        # Carbons can be easily determined by their # bonds
        atm = molecule.atms[iatm]
        nbonds = len(atm.bonds)
        if atm.hyb == 9:
            prefix = 'CR'
        elif nbonds == 4:
            prefix = 'CS'
        elif nbonds == 3:
            prefix = 'CD'
        elif nbonds == 2:
            prefix = 'CT'
        else:
            return False

        nH = 0
        for jatm,order in atm.bonds:
            if molecule.atms[jatm].is_H:
                nH += 1

        # make this decision afterwards...
        #if atm.connected_to_polar:
        #    if nH < 2: # keep original if more than two hydrogens attached
        #        return prefix + 'p'
        #    else:
        #        print "Assign non-polar based on nH: ", atm.name, prefix
            
        if prefix in ['CS','CD']:
            atype = prefix
            if nH > 0:
                atype += '%d'%nH
            if prefix not in ACLASS_ID: 
                atype = 'Null'
            return atype
        return prefix

    def classify_O(self,molecule,iatm):
        atm = molecule.atms[iatm]

        # first check if is connected to any carbon
        atoms_connected = [molecule.atms[jatm].atype for jatm,order in atm.bonds]
        nH = atoms_connected.count(2)

        # return generic if not connected to C
        if atoms_connected.count(1) == 0:
            if atm.hyb == 3:
                if nH >= 1: #but still should be 1
                    is_PO4H = False
                    if len(atm.bonds) == 2:
                        for jatm,order in atm.bonds:
                            if molecule.atms[jatm].atype == 5:
                                is_PO4H = True
                                break
                    if is_PO4H:
                        return 'Ohx' #PO4H: should revisit someday
                    else:
                        return 'OG31' #mostly for OXIME now
                else:
                    return 'OG3'
            elif atm.hyb == 2:
                # nitro group
                if atoms_connected.count(4) == 1 and len(atm.bonds) == 1:
                    #print 'detected nitro'
                    return 'Ont'
                else:
                    return 'OG2'
                
            #elif atm.hyb == 'R':
            # could this happen?    
            else: #shouldn't happen
                #print 'Warning: not connected to carbon but hybridization neither 2 nor 3'
                return False

        # from here should be connected to at least one C
        # check if nbonds <=2
        if len(atm.bonds) > 2:
            return False

        # anything aromatic: simplest -- this doesn't occur unfortunately; there is no O.ar
        #if atm.hyb == 9:
        #    return 'Ofu'

        # sp3
        elif atm.hyb == 3:
            #hydroxyl, carboxylic acid: anything C-O-H
            if atoms_connected.count(2) == 1: 
                return 'Ohx'
            #ether, ester, acetyl, furan
            elif len(atm.bonds) == 2 and atoms_connected.count(1) == 2: #cannot be >2
                #nCaro = 0
                #for jatm,order in atm.bonds:
                #    if molecule.atms[jatm].hyb in [2,9]:
                #        nCaro += 1
                #if nCaro == 2: # as 5 or 6-membered ring
                # modified Jan 2019
                if iatm in molecule.atms_aro:
                    return 'Ofu'
                else:
                    return 'Oet'
            else:  #generic
                if nH >= 1:
                    return 'OG31'
                else:
                    return 'OG3'

        # sp2
        # should look to special case of amide
        elif atm.hyb == 2: 
            if atoms_connected.count(1) == 2:
                # this could only happen for Furan, when hyb of O is set to 2
                return 'Ofu'

            atm_j = molecule.atms[ atm.bonds[0][0] ] # has only 1 bond!
            # make sure this is carbon; cannot be other 
            if atm_j.atype != 1: 
                #print 'Warning: atm_j is not carbon but ', atm_j.atype
                #return False
                return 'OG2'

            atoms_connected_to_j = [molecule.atms[katm].atype for katm,order in atm_j.bonds]
            hybs_k = [molecule.atms[katm].hyb for katm,order in atm_j.bonds]
            nC = atoms_connected_to_j.count(1)
            nH = atoms_connected_to_j.count(2)
            nO = atoms_connected_to_j.count(3)-1 # exclude self
            nN = atoms_connected_to_j.count(4)

            if nN >= 1:
                return 'Oad' # amide; allow N-C(O)-N as amide
            #elif nN == 2:
            #    return 'OG2' # O attached to -N-C-N-; is this same as amide or not? (see up)

            elif nN == 0: 
                if nC  == 2:
                    return 'Oal' #ketone
                elif nC + nH == 2: # nC/nH = 1/1 or 0/2, aldehyde
                    return 'Oal'

                #   O
                # R-C-O
                elif nO == 1:
                    # first get the other O
                    for bond in atm_j.bonds:
                        if bond[0] == iatm: continue #skip itself
                        atm_k = molecule.atms[bond[0]]
                        if atm_k.atype == 3:
                            break
       
                    # look how many bonds the other O has: -O or -O-X
                    nbonds_k = len(atm_k.bonds)
                    if nbonds_k == 2: 
                        return 'Oal' # O=C-O-X: acetyl, carboxylate, ester...
                    elif nbonds_k == 1:
                        return 'Oat' # O=C-O
                    else:
                        return 'OG2'
                else:
                    return 'OG2'
            else:
                return 'OG2'
                #return False

    def classify_N(self,molecule,iatm):
        atm = molecule.atms[iatm]

        # special rule for sp(1)
        if atm.hyb == 1:
            return 'NG1'

        # first check if is connected to any carbon
        atoms_connected = [molecule.atms[jatm].atype for jatm,order in atm.bonds]
        nC = atoms_connected.count(1)
        nH = atoms_connected.count(2)
        nO = atoms_connected.count(3)
        nN = atoms_connected.count(4)
        ntot = len(atm.bonds)

        # 1. Consider case if something else (including lone pair(s)) than C/H connected 
        # exceptionally allow amide to pass here even if O=C-N-X where X = C or H
        if nC+nH < ntot:
            if atm.hyb == 3:
                if len(atm.bonds) <= 3 and nH >= 1: 
                    return 'Nam2' # Amine with lone pair(s)
                elif nH >= 1:
                    return 'Nam'
                else: 
                    return 'NG3' # sp3, 0H
                
            elif atm.hyb in [2,8,9]: #sp2 or aromatic
                if nO >= 2:
                    if nH >= 1:
                        return False
                    return 'NG2' # nitro (share atomtype; oxygen is more important)

                #elif nS >= 1: #Sulfonate: currently NG2, but seperate it?
                elif nH == 0 and nN >= 1: # Azo group; -N=N- 
                    return 'Nad3'
                else:
                    if nH == 0:
                        return 'NG2'
                    elif nH == 1:
                        #explicitly rescue possible amides...
                        is_amide = False # is_amide = (atm.hyb==8) 
                        for jatm,jorder in atm.bonds: 
                            atm_j = molecule.atms[jatm]
                            if atm_j.atype == 1 and atm_j.hyb == 2: #not aromatic
                                for katm,korder in atm_j.bonds: 
                                    atm_k = molecule.atms[katm]
                                    if atm_k.atype == 3 and atm_k.hyb == 2:
                                        is_amide = True
                                        break
                        if is_amide:
                            return 'Nad'
                        else:
                            return 'NG21'
                    else:
                        return 'NG22'
                    
            else: #shouldn't happen
                return False

        # From here should be connected to H or C only, no lone pair
        # Amide:
        # Caution: sometimes hyb can be else than 8 depending on mol2 preparation
        elif atm.hyb == 8: #amide
            if nC == 3:
                return 'Nad3'
            else: # primary/secondary amide
                if nH == 0:
                    return 'Nad3'
                else:
                    return 'Nad'
            
        # Anything aromatic
        elif atm.hyb == 9:
            if ntot == 3 and nC == 2 and nH == 1:
                return 'Nin' # Iminium or Indole
            elif ntot == 3 and nC == 3: 
                return 'Nad3' #C=NC=C
            elif ntot == 2 and nC == 2:
                return 'Nim' # Imidazole; deprotonated
            else:# O=N-Caro, Naro-Naro-Caro, ...
                if nH == 0:
                    return 'NG2' 
                elif nH == 1:
                    return 'NG21'
                else:
                    return 'NG22'

        # sp3
        elif atm.hyb == 3:
            if nH == 0:
                #if len(atm.bonds) <= 3: #treat lone-pair case separately?
                return 'NG3' 
            else:
                if len(atm.bonds) <= 3:
                    return 'Nam2' #amine/ammonium with a lone-pair
                else: 
                    return 'Nam' #amine/ammonium

        elif atm.hyb == 2:
            # guadinium / amide precedes anything else if condition satisfied
            is_guanidinium = False
            is_amide = False
            nCaro = 0 

            # Iter through connected atoms
            for jatm,jorder in atm.bonds: 
                atm_j = molecule.atms[jatm]
                # Count no. Caro or Csp2 
                if atm_j.atype == 1 and atm_j.hyb in [2,9]:
                    nCaro += 1

                # Functional group with sp2 carbon: guanidinium
                if (atm_j.atype == 1 and atm_j.hyb in [2,8,9]):
                    nOsp2_j = 0 #Osp2 at 1-3; consider as amide
                    nNsp2_j = 0 #Nsp2 at 1-3; consider as guanidinium if any j attached to 2 of these
                    #nCp_j = 0
                    
                    # iter through 1-3 pairs
                    for katm,korder in atm_j.bonds:
                        atm_k = molecule.atms[katm]
                        if katm == iatm: continue

                        if atm_k.atype == 3 and atm_k.hyb == 2: #Osp2
                            nOsp2_j += 1
                        elif atm_k.atype == 4 and atm_k.hyb in [2,8,9]: # Nsp2
                            nNsp2_j += 1
                        #elif atm_k.atype == 1 and ACLASS_ID[atm.aclass] in ['CDp','CRp']: # nCpolar
                        #    nCp_j += 1

                    if nOsp2_j == 1:
                        is_amide = True
                    if nNsp2_j == 2: # if any carbon has two sp2 nitrogen attached
                        is_guanidinium = True

            #print 'Here!', is_amide, nC, nOsp2
            if is_amide: #overrides else
                if nC == 3: # tertiary, 0-H 
                    return 'Nad3'
                elif nH >= 1: # primary/sec, donor
                    return 'Nad'
                else: # primary/secondary amide, non-donor
                    return 'NG2' # share with Naro_0H; similar property?

            if is_guanidinium:
                if nH == 2: 
                    return 'Ngu2' #R-NH2
                elif nH == 1:
                    return 'Ngu1'# R-NH-Gu
                else:
                    # R-NR-Gu or R-N:-Gu; share with Naro_0H; similary property?
                    # TODO: try Nim for acceptor property
                    return 'NG2' 
                
            elif nC == 2 and nCaro in [1,2]: #attached to at least one Caro
                if nH == 1 and iatm in molecule.atms_aro:
                    return 'Nin' #Oct24/2018 Nin is only for non-aromatic
                elif nH == 1:
                    return 'NG21'
                elif nH == 0:
                    return 'Nim'
                
            elif nC == 3 and nCaro >= 1: # planary C=NC=C, 0H
                return 'Nad3' # assign as tertiary amide
                
            else:
                if nH == 0:
                    return 'NG2' 
                elif nH == 1:
                    return 'NG21'
                else: #sp2 cannot have > 2H
                    return 'NG22'
        else: 
            return False

    def classify_S(self,molecule,iatm):
        atm = molecule.atms[iatm]

        # first check if is connected to any carbon
        atoms_connected = [molecule.atms[jatm].atype for jatm,order in atm.bonds]
        nC = atoms_connected.count(1)
        nH = atoms_connected.count(2)
        nS = atoms_connected.count(6)
        ntot = len(atoms_connected)
        if nC == 1 and nH == 1 and ntot == 2:
            return 'Sth' #thiol
        #elif nC == 1 and ntot == 1: # this is unlikely, but because of artifact in CSD?
        #    return 'Sth' #thiol
        elif nC+nS == 2 and ntot == 2:
            if iatm in molecule.atms_aro:
                return 'SR' #aromatic; Jan 2019
            else:
                return 'Ssl' #sulfide/disulfide
        elif ntot == 1:
            #print 'detected SG2'
            return 'SG2'
        else:
            if atm.hyb == 5:
                return 'SG5'
            #elif atm.hyb == 3:
            #    return 'SG3'
            #elif atm.hyb == 2:
            #    return 'SG2'
            #else:
            #    return False
            else:
                return 'SG3'

    def classify_P(self,molecule,iatm):
        atm = molecule.atms[iatm]
        if atm.hyb == 5:
            return 'PG5'
        else:
            return 'PG3'

    def classify_halogen(self,molecule,iatm):
        atm = molecule.atms[iatm]
        prefix = ATYPES_REG[atm.atype]
        if len(atm.bonds) == 1:
            j = atm.bonds[0][0]
            atm_j = molecule.atms[ j ]
            #if atm_j.atype == 1 and atm_j.hyb == 9:
            if atm_j.atype == 1 and j in molecule.atms_aro:
                prefix += 'R'
        return prefix

    def classify_else(self,molecule,iatm):
        # special atms
        if molecule.atms[iatm].atype >= len(ATYPES_REG)-1:
            aname = ATYPES[molecule.atms[iatm].atype]
            if aname in ['Ca','Mg','Fe','Zn','Co','Cu']:
                return aname+'2p'
            elif aname in ['B']:
                return aname+ 'sp2'
            else: #Si,Se,Mn,Cd,Ni
                return aname
        else:
            print('Unknow atom: ', molecule.atms[iatm].report_self())
            return False

    def modify_polarC(self,molecule,iatm):
        atm = molecule.atms[iatm]

        if atm.atype != 1: return
        aclass = atm.aclass

        nH = 0
        for jatm,order in atm.bonds:
            if molecule.atms[jatm].is_H:
                nH += 1
        
        if len(atm.bonds) - nH == 1: # pass if is at aliphatic tip
            return
            
        attached_to_polar = False
        for jatm,order in atm.bonds:
            jclass = ACLASS_ID[molecule.atms[jatm].aclass]
            if jclass in POLARCLASSES:
                attached_to_polar = True
                break
            
        if attached_to_polar:
            molecule.atms[iatm].aclass = ACLASS_ID.index(ACLASS_ID[aclass][:2]+'p')
            #print 'Modified atype to polar: ', atm.name, aclass, molecule.atms[iatm].aclass
        #if (not attached_to_polar) and abs(atm.charge) > 0.3:
        #    print('Warning: not assigned as polar carbon but still has some partial charge: ', atm.name, ACLASS_ID[aclass], atm.charge)

    def assert_H(self,molecule):
        errormsg = ''
        for iatm,atm in enumerate(molecule.atms):
            aclass = ACLASS_ID[atm.aclass]
            if aclass not in ASSERT_NHYDROGENS: continue

            nH = 0
            for jatm,order in atm.bonds:
                if molecule.atms[jatm].is_H:
                    nH += 1
            nH_assert = ASSERT_NHYDROGENS[aclass]
            if nH not in nH_assert:
                errormsg += '%2s %4s %1d != %s\n'%(atm.name, aclass, nH, ','.join(['%d'%n for n in nH_assert]))

        if errormsg != '':
            sys.exit('Error, inconsistent num hydrogens!\n'+errormsg)
        return

class FunctionalGroupClassifier:
    def __init__(self):
        return

    def apply_to_molecule(self,molecule):
        WCCLASS = {'C*':['CS','CS1','CS2','CS3','CSp','CD','CD1','CD2','CDp',
                         'CR','CRp','CT','CTp'],
                   'CS*':['CS','CS1','CS2','CS3','CSp'],
                   'CD*':['CD','CD1','CD2','CDp'],
                   'CR*':['CR','CRb','CRp'],
                   'OG2*':['Ohx','Oet','OG3','OG31'],
                   'OG3*':['Oad','Oat','Ofu','Ont','OG2'],
                   'O*':['Ohx','Oet','Oal','Oad','Oat','Ofu','Ont','OG2','OG3','OG31'],
                   'NG2*':['Nam','Nam2','NG3'],
                   'NG3*':['Nad3','Nin','Nim','Ngu1','Ngu2','NG2', 'NG21','NG22','NG1'],
                   'N*':['Nam','Nam2','Nad','Nad3','Nin','Nim','Ngu1','Ngu2','NG3','NG2', 'NG21','NG22','NG1'],
                   'X':['F','Br','Cl','I'],
                   'XR':['FR','BrR','ClR','IR']
                   }


        TORSKEY2FUNCGRP = \
        {#'CR.CRb.CRb.CR':'biaryl',
         #'CR.CRb.CRb.CD':'biaryl',
         #'CR.CRb.NGb.CR':'biaryl',
         #'CD.CRb.NGb.CR':'biaryl',
         #'CD.CRb.NGb.CD':'biaryl',
            'CS*.CS2.CS2.CS*':'aliphatic',
            'HO.Ohx.CDp.Oal':'carboxylicacid',
            #'Oet.CS2.Oet.CS*':'acetal', # can precede ether
            #'Oet.CS1.Oet.CS*':'acetal', # can precede ether 
            #'Oet.CS.Oet.CS*':'acetal', # can precede ether
         'CR.Oet.CDp.Oal':'aryl-ester',
         'CS*.Oet.CDp.Oal':'ester',
         'CD.CD.Oet.CS*':'aryl-ester', #precedes regular
         'CR.CR.Oet.CS*':'aryl-ether', #precedes other aryls
         'CR.CRp.Ohx.HO':'aryl-hydroxyl', #precedes other aryls
         'CR.CR.NG2.Ont':'aryl-nitro', #precedes regular nitro
         'CR.CR.CDp.Oad':'aryl-amide',
         'CR.CR.Nad.CDp':'aryl-amide',
         'CD*.CDp.NG22.HN':'aniline',
         #'C*.CRp.NG21.HN':'aryl-2amine',
            #'CD.CDp.Nin.HN':'aryl-2amine',
            #'C*.CRp.Nin.HN':'aryl-2amine',
            #'C*.CR.Nin.HN':'aryl-2amine',
         'Oad.CDp.Nad.HN':'amide',
         'Oad.CDp.Nad3.CS*':'amide3',
         'C*.NG2.OG31.HO':'oxime',
         'C*.Nad3.Nad3.C*':'azo',
         'C*.Nad.Nad3.C*':'azo',
         'C*.NG21.Nad3.C*':'azo',
         'CR.CRp.SG5.OG2':'aryl-sulfonate', #precedes regular
         'CR.Nam2.SG5.OG2':'aryl-sulfonamide', #precedes regular
         'CR.NG21.SG5.OG2':'aryl-sulfonamide', #precedes regular
        }
        
        ANGKEY2FUNCGRP = \
        {'CS*.Oet.CS*':'ether',
         'CD*.Oet.CS*':'ether',
         'Oat.CDp.Oat':'carboxylate',
         'CR.Ofu.CR':'furan',
         'CD*.Ofu.CD*':'furan',
         'Oal.CDp.CS*':'ketone/aldehyde',
         'Oal.CD.CS*':'ketone/aldehyde',
         'CS*.Ohx.HO':'alcohol',
         'Ont.NG2.Ont':'nitro',
         'C*.Nam.HN':'amine',
         'C*.Nam2.HN':'amine3',
         'HN.NG22.HN':'amine',
         #'C*.CT.NG1':'nitrile',
         'C*.Sth.HS':'thiol',
         'C*.CD.SG2':'sulfone', #check
         'OG2.SG5.OG2':'sulfonate',
         'OG2.PG5.OG2':'phosphate',
         'C*.CS.X':'aliphatichalide',
         'C*.CS1.X':'aliphatichalide',
         'C*.CS2.X':'aliphatichalide',
         'C*.CS3.X':'aliphatichalide',
         'CR.CR.XR':'aryl-halide',
         'CR.CR.NG2*':'aryl',
         'CR.CR.CD1':'aryl',
         'CR.CR.CDp':'aryl',
         'SG2.CDp.N*':'thio-', #include -amide
         'SG2.CDp.C*':'thio-', #include -ketone..
         'SG2.CD.N*':'thio-', #include -amide
         'SG2.CD.C*':'thio-', #include -ketone..
         'HN.Ngu2.CDp':'guanidinium',
         'CR.Nim.CR':'pyridine',
         'C*.NG1.NG1':'diazo',
         }

        BNDKEY2FUNCGRP = \
        {'CT.NG1':'nitrile'
        }

        GRPS_NOT_IN_RING = ['ether','ester','aryl',]
        ALLOW_MULTI_DEF = ['amide','azo']#,'aryl-sulfonate','aryl-sulfonamide']
        funcgrp_assigned = []

        # biaryl pivots
        for pivot in molecule.biaryl_pivots:
            atms = tuple([molecule.atms[i] for i in pivot])
            if pivot in molecule.biaryl_pivots_extra:
                funcgrp = FunctionalGroupClass(atms,"BiarylExtra")
                molecule.funcgrps.append(funcgrp)
            else:
                funcgrp = FunctionalGroupClass(atms,"BiarylRegular")
                molecule.funcgrps.append(funcgrp)
        
        # special logic for sugar or imidazole/pyrrole/indole
        for ring in molecule.rings:
            #if ring in molecule.rings_aro: continue
            if ring.natms > 6 or ring.natms < 5: continue

            has_Oet = False
            has_Ohx = False
            has_nonsp3_carbon = False
            has_sp3_carbon = False
            #has_Nim = False
            has_Nin = False
            has_NN = False
            has_Ssl = False
            for a in ring.atms:
                t = ACLASS_ID[molecule.atms[a].aclass]
                is_nonsp3_carbon = False
                for key in ['CD*','CR*']:
                    if t in WCCLASS[key]:
                        has_nonsp3_carbon = True
                        break
                if t in WCCLASS['CS*']:
                    has_sp3_carbon = True
                if t in ['Oet','Ofu']:
                    has_Oet = True
                elif t == 'Nin':
                    has_Nin = True
                elif t == 'Ssl':
                    has_Ssl = True
                else:
                    for b,order in molecule.atms[a].bonds:
                        if molecule.atms[b].is_H: continue
                        if b in ring.atms:
                            if t == 'Nad3' and ACLASS_ID[molecule.atms[b].aclass] == 'Nad3':
                                has_NN = True
                        else:
                            tb = ACLASS_ID[molecule.atms[b].aclass]
                            if tb == 'Ohx':
                                has_Ohx += 1
            atms = tuple([molecule.atms[i] for i in ring.atms])
            if has_NN:
                funcgrp = FunctionalGroupClass(atms,'pyrazole')
                molecule.funcgrps.append(funcgrp)
                funcgrp_assigned += ring.atms
            elif has_Oet and has_Ohx and not has_nonsp3_carbon and (ring.type != 3):
                funcgrp = FunctionalGroupClass(atms,'sugar')
                molecule.funcgrps.append(funcgrp)
                funcgrp_assigned += ring.atms
            elif has_Nin and (ring.type != 3):
                funcgrp = FunctionalGroupClass(atms,'pyrrole')
                molecule.funcgrps.append(funcgrp)
                funcgrp_assigned += ring.atms
            elif has_Ssl and (ring.type != 3) and not has_sp3_carbon:
                funcgrp = FunctionalGroupClass(atms,'thiolane')
                molecule.funcgrps.append(funcgrp)
                funcgrp_assigned += ring.atms
                        
        # first based on torsions
        for i,torsion in enumerate(molecule.torsions):
            assigned = False
            in_ring = 0
            for a in torsion:
                if a in funcgrp_assigned:
                    assigned = True
                if a in molecule.atms_ring:
                    in_ring += 1
            if assigned: continue
            in_ring = (in_ring > 2)
            (i1,i2,i3,i4) = torsion

            t1 = ACLASS_ID[molecule.atms[i1].aclass]
            t2 = ACLASS_ID[molecule.atms[i2].aclass]
            t3 = ACLASS_ID[molecule.atms[i3].aclass]
            t4 = ACLASS_ID[molecule.atms[i4].aclass]

            key = '%s.%s.%s.%s'%(t1,t2,t3,t4)
            key_inv = '%s.%s.%s.%s'%(t4,t3,t2,t1)
            keys = [key,key_inv]

            wc1s = []
            wc4s = []
            for wc in WCCLASS:
                if t1 in WCCLASS[wc]:
                    wc1s.append(wc)
                if t4 in WCCLASS[wc]:
                    wc4s.append(wc)

            for wc1 in wc1s:
                keys.append('%s.%s.%s.%s'%(wc1,t2,t3,t4))
                keys.append('%s.%s.%s.%s'%(t4,t3,t2,wc1))
                for wc4 in wc4s:
                    keys.append('%s.%s.%s.%s'%(wc1,t2,t3,wc4))
                    keys.append('%s.%s.%s.%s'%(wc1,t3,t2,wc4))
            for wc4 in wc4s:
                keys.append('%s.%s.%s.%s'%(t1,t2,t3,wc4))
                keys.append('%s.%s.%s.%s'%(wc4,t3,t2,t1))

            #print (key, keys)
            for key in keys:
                if key in TORSKEY2FUNCGRP:
                    grpname = TORSKEY2FUNCGRP[key]
                    if in_ring and (grpname in GRPS_NOT_IN_RING): continue
                    funcgrp = FunctionalGroupClass((molecule.atms[i1],molecule.atms[i2],molecule.atms[i3],molecule.atms[i4]),
                                                   grpname)
                    molecule.funcgrps.append(funcgrp)
                    #if grpname in ALLOW_MULTI_DEF: continue
                    funcgrp_assigned += [i2,i3]
                    #else:
                    #funcgrp_assigned += [i1,i2,i3,i4]
                    break
            
        # then based on angle
        for i,angle in enumerate(molecule.angles):
            assigned = False
            in_ring = 0
            for a in angle:
                if a in funcgrp_assigned: assigned = True
                if a in molecule.atms_ring:
                    in_ring += 1
            if assigned: continue
            in_ring = (in_ring > 1)
            (i1,i2,i3) = angle
            
            t1 = ACLASS_ID[molecule.atms[i1].aclass]
            t2 = ACLASS_ID[molecule.atms[i2].aclass]
            t3 = ACLASS_ID[molecule.atms[i3].aclass]

            key = '%s.%s.%s'%(t1,t2,t3)
            key_inv = '%s.%s.%s'%(t3,t2,t1)
            keys = [key,key_inv]

            wc1s = []
            wc2s = []
            wc3s = []
            for wc in WCCLASS:
                if t1 in WCCLASS[wc]:
                    wc1s.append(wc)
                if t2 in WCCLASS[wc]:
                    wc2s.append(wc)
                if t3 in WCCLASS[wc]:
                    wc3s.append(wc)

            #print (molecule.atms[i1].name,t1, wc1s)
            #print (molecule.atms[i3].name,t3, wc3s)
            for wc1 in wc1s:
                keys.append('%s.%s.%s'%(wc1,t2,t3))
                keys.append('%s.%s.%s'%(t3,t2,wc1))
                for wc3 in wc3s:
                    keys.append('%s.%s.%s'%(wc1,t2,wc3))
                    keys.append('%s.%s.%s'%(wc3,t2,wc1))
            for wc3 in wc3s:
                keys.append('%s.%s.%s'%(t1,t2,wc3))
                keys.append('%s.%s.%s'%(wc3,t2,t1))
            #print(keys)
            
            for key in keys:
                if key in ANGKEY2FUNCGRP:
                    grpname = ANGKEY2FUNCGRP[key]
                    if in_ring and (grpname in GRPS_NOT_IN_RING): continue
                    funcgrp = FunctionalGroupClass((molecule.atms[i1],molecule.atms[i2],molecule.atms[i3]),
                                                   grpname)
                    molecule.funcgrps.append(funcgrp)
                    funcgrp_assigned += [i1,i2,i3]
                    break

        # finally based on angle
        for i,bond in enumerate(molecule.bonds):
            i1 = bond.atm1
            i2 = bond.atm2
            if i1 in funcgrp_assigned: continue
            if i2 in funcgrp_assigned: continue
            
            t1 = ACLASS_ID[molecule.atms[i1].aclass]
            t2 = ACLASS_ID[molecule.atms[i2].aclass]

            key = '%s.%s'%(t1,t2)
            key_inv = '%s.%s'%(t2,t1)
            keys = [key,key_inv]

            wc1s = []
            wc2s = []
            for wc in WCCLASS:
                if t1 in WCCLASS[wc]:
                    wc1s.append(wc)
                if t2 in WCCLASS[wc]:
                    wc2s.append(wc)

            for wc1 in wc1s:
                keys.append('%s.%s'%(wc1,t2))
                keys.append('%s.%s'%(t2,wc1))
                for wc2 in wc2s:
                    keys.append('%s.%s'%(wc1,t2))
                    keys.append('%s.%s'%(t1,wc2))
            for wc2 in wc2s:
                keys.append('%s.%s'%(t1,wc2))
                keys.append('%s.%s'%(wc2,t1))

            for key in keys:
                if key in BNDKEY2FUNCGRP:
                    funcgrp = FunctionalGroupClass((molecule.atms[i1],molecule.atms[i2]),
                                        BNDKEY2FUNCGRP[key])
                    molecule.funcgrps.append(funcgrp)
                    funcgrp_assigned += [i1,i2]
                    break

    
