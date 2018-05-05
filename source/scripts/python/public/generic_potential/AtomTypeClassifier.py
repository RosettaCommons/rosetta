import sys

from Types import *

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
                    return 'OG31'
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
                nCaro = 0
                for jatm,order in atm.bonds:
                    if molecule.atms[jatm].hyb in [2,9]:
                        nCaro += 1
                if nCaro == 2: # as 5 or 6-membered ring
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
                if nH == 1:
                    return 'Nin'
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
            atm_j = molecule.atms[ atm.bonds[0][0] ]
            if atm_j.atype == 1 and atm_j.hyb == 9:
                prefix += 'R'
        return prefix

    def classify_else(self,molecule,iatm):
        if molecule.atms[iatm].atype >= len(ATYPES_REG)-1:
            aname = ATYPES[molecule.atms[iatm].atype]
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
        if (not attached_to_polar) and abs(atm.charge) > 0.3:
            print('Warning: not assigned as polar carbon but still has some partial charge: ', atm.name, ACLASS_ID[aclass], atm.charge)

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
