import sys
import operator
from utils import distance,angle,dihedral
from Types import *

class AtomClass:
    def __init__(self,name,atype,hyb,charge):
        self.hyb = hyb
        self.name = name
        self.atype = atype
        self.bonds = []
        self.connected_to_polar = False
        self.has_H = False
        self.is_H = (atype == 2)
        self.aclass = 0
        self.charge = charge
        self.icoord = []
        self.root = -1 # undefined
    
    def add_bond(self,atm_connected,order):
        self.bonds.append((atm_connected,order))

    def report_self(self):
        return ' %2s %6s %3d'%(self.atype,self.name,self.hyb)

class BondClass:
    def __init__(self,atm1,atm2,order):
        self.atm1 = min(atm1,atm2)
        self.atm2 = max(atm1,atm2)
        self.order = order

class MoleculeClass:
    def __init__(self,mol2file,verbose=False,):
        self.name = ''
        self.atms = []
        self.bonds = []
        self.angles = []
        self.torsions = []
        self.nbratom = 0
        self.nbrradius = 0.0
        self.torsionparams = []
        self.torsion_orders = []
        self.hapol_torsion_id = []
        self.hpol_torsion_type = []
        self.rotable_torsion_id = []
        self.xyz = []
        self.nheavyatm = 0
        self.atms_aro = []
        self.atms_puckering = [] #including aromatic ring
        self.max_confs = 5000 #copied from default of molfile_to_params
        self.crystinfo = ''
        self.rings_aro = [] #only aromatic
        self.rings = [] #aromatic + puckering
        self.chargegrps = []

        self.chiatms = []
        self.chitypes = []
        self.chiextra = ''
        
        self.mol2file = mol2file
        stat = self.read_mol2(verbose=verbose)
        if not stat:
            return
        self.assign_angles()
        self.assign_torsions(verbose=verbose)

    def read_mol2(self,verbose=False):
        mode = 0
        i_mode1 = 0
        newindex = []
        atms = []
        xyzs = []
        bonds = []
        #self.atomnames_in_inputmol2 = []

        for l in open(self.mol2file):
            words = l[:-1].split()
            if l.startswith('@<TRIPOS>MOLECULE'):
                mode = 1
                continue
            elif l.startswith('@<TRIPOS>ATOM'):
                mode = 2
                continue
            elif l.startswith('@<TRIPOS>BOND'):
                mode = 3
                continue
            elif l.startswith('@<TRIPOS>CRYSIN'):
                mode = 4
                continue

            if mode == 1:
                i_mode1 += 1
                if i_mode1 == 1:
                    self.name = words[0]

            if mode == 2 and len(words) == 9:
                i = int(words[0])
                atype1 = words[5].split('.')[0].upper()
                #if atype in ATYPES:
                #    atype = atype2

                if atype1 in ATYPES:
                    pass
                else:
                    atype1 = self.try_map_atype(atype1)
                    
                atype = ATYPES.index(atype1) # integer
                name = words[1] # string

                if '?' in name: # ambiguous
                    continue
                    
                hybstr = words[5]
                hyb = 0 # integer

                if hybstr in SPECIAL_HYBRIDS:
                    hyb = SPECIAL_HYBRIDS[hybstr]
                elif '.' in hybstr:
                    hyb_i = hybstr.split('.')[-1]
                    if hyb_i in ['1','2','3']:
                        hyb = int(hyb_i)
                    elif '.ar' in hybstr:
                        hyb = 9
                    else:
                        hyb = 0
                else:
                    hyb = 0 # none

                charge = float(words[-1])

                xyz = [float(word) for word in words[2:5]]
                atms.append( AtomClass(name,atype,hyb,charge) )
                xyzs.append( xyz )
                newindex.append((atype==2,len(atms)-1))
                #self.atomnames_in_inputmol2.append(name)
                #atms_sortable.append( (atype==2,AtomClass(name,atype,hyb,charge),xyz) )

            elif mode == 3 and len(words) == 4:
                i = int(words[0])
                atm1 = int(words[1])-1 #0-index
                atm2 = int(words[2])-1 #0-index
                if words[3] in ['1','2','3']:
                    order = int(words[3])
                elif words[3] == 'am':
                    order = 2 #amide
                elif words[3] == 'ar':
                    order = 9 #aromatic
                #elif words[3] == 'un': 
                #    #order = 0 #
                #    print '
                #  
                else:
                    print('skip bond order of ', words[3])
                    order = 0
                    
                bonds.append((atm1,atm2,order))
            elif mode == 4:
                self.parse_crystinfo(l)
                mode = 0

        # Re-ordering scheme for Rosetta
        # add up here with new index based on heavy or H
        newindex.sort()
        newindex = [comp[1] for comp in newindex]
        self.atms = []
        self.xyz = []
        for inew,iprv in enumerate(newindex):
            if verbose:
                print('Sorted atom index: %s %d -> %d'%(atms[iprv].name,iprv,inew))
            self.atms.append(atms[iprv])
            self.xyz.append(xyzs[iprv])
            if not self.atms[inew].is_H: self.nheavyatm = inew+1

        bonds_reindexed = []
        for atm1,atm2,order in bonds:
            atm1n = newindex.index(atm1)
            atm2n = newindex.index(atm2)
            bonds_reindexed.append((
                (self.atms[atm1n].is_H or self.atms[atm2n].is_H),
                atm1n,atm2n,order))

        bonds_reindexed.sort()

        self.bonds = []
        for has_H,atm1n,atm2n,order in bonds_reindexed:
            self.bonds.append(BondClass(atm1n,atm2n,order))
        self.bonds.sort(key=operator.attrgetter('order'))

        # assign bonds
        self.assign_bonds()

        # assign hybridization in case not specified in mol2
        self.assign_hybridization()
        return True

    # I don't recommend this... use already assigned one
    def assign_hybridization(self):
        n_to_assign = 0
        to_assign = []
        for iatm,atm in enumerate(self.atms):
            if atm.atype in ATYPES_HYBRID:
                n_to_assign += 1
                if atm.hyb == 0:
                    print('WARNING: Hybridization state for %d:%s is not assigned!'%(iatm,atm.name))
                    to_assign.append(iatm)

        return
        for iatm in to_assign:
            atm = self.atms[iatm]
            nbonds = sum([bond.bondorder for bond in atm.bonds])
            if atm.atype == 1: #Carbon
                if nbonds == 4: self.atms[iatm].hyb = 3
                elif nbonds == 3: self.atms[iatm].hyb = 2
                elif nbonds == 2: self.atms[iatm].hyb = 1
            elif atm.atype == 3: #Oxygen
                if nbonds == 2: self.atms[iatm].hyb = 3
                elif nbonds == 1: self.atms[iatm].hyb = 2
            elif atm.atype == 4: #Nitrogen
                if nbonds == 4: self.atms[iatm].hyb = 3
                elif nbonds == 3: self.atms[iatm].hyb = 2
                elif nbonds == 2: self.atms[iatm].hyb = 1
            
    def try_map_atype(self,atype):
        if atype.upper() == 'BR':
            return 'Br'
        elif atype.upper() == 'CL':
            return 'Cl'
        elif atype[0] in ATYPES_REG:
            return atype[0]
        else:
            print('No such atom type: %s, while reading mol2file '%atype+self.mol2file+', line: ')
            print()
            print(l)
            return 'Null'

    def assign_bonds(self):
        # first, correction on mol2 file...
        #for i,bond in enumerate(self.bonds):
        #    if self.atms[bond.atm1].hyb > 1 and self.atms[bond.atm2].hyb > 1:
        #        if bond.order == 1: 
        #            self.bonds[i].order = 2

        # then go through misc...
        for bond in self.bonds:
            is_H1 = (self.atms[bond.atm1].atype==2) 
            is_H2 = (self.atms[bond.atm2].atype==2) 
            self.atms[bond.atm1].add_bond(bond.atm2,bond.order)
            self.atms[bond.atm2].add_bond(bond.atm1,bond.order)

            if is_H1:
                self.atms[bond.atm2].has_H = True
            if is_H2:
                self.atms[bond.atm1].has_H = True

            if self.atms[bond.atm1].atype == 1 and self.atms[bond.atm2].atype in POLAR_ATOMS:
                self.atms[bond.atm1].connected_to_polar = True
            if self.atms[bond.atm2].atype == 1 and self.atms[bond.atm1].atype in POLAR_ATOMS:
                self.atms[bond.atm2].connected_to_polar = True

        # finally, sanity check if any atm is not connected to anything
        atms_unconnected = []
        for atm in self.atms:
            if len(atm.bonds) == 0:
                atms_unconnected.append(atm.name)
        if atms_unconnected != []:
            sys.exit('ERROR: Atom found not connected to any other atom:'+' '.join(atms_unconnected))
            


    def assign_angles(self):
        for i,bond1 in enumerate(self.bonds[:-1]):
            for bond2 in self.bonds[i+1:]:
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
                self.angles.append(angle)

    def assign_torsions(self,verbose=False):
        bonds = [(bond.atm1,bond.atm2) for bond in self.bonds]

        torsions_heavy = []
        torsions_H = []
        orders_heavy = []
        orders_H = []
        if verbose:
            print('bonds: ', [(bond.atm1,bond.atm2) for bond in self.bonds])
        for i,bond1 in enumerate(self.bonds[:-1]):
            orders_tmp = []
            torsions_tmp = []
            atypes = []
            dtypes = []
            for bond2 in self.bonds[i+1:]:
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
                orders_tmp.append(self.bonds[ibond].order)
                torsions_tmp.append((a,b,c,d))
                
                atypes.append(self.atms[a].is_H)
                dtypes.append(self.atms[d].is_H)

            #if atypes.count(True) == len(atypes) or dtypes.count(True) == len(dtypes): #if having only hydrogens at either end
            if atypes.count(True) > 0 or dtypes.count(True) > 0: #if having any hydrogens at either end
                orders_H += orders_tmp
                torsions_H += torsions_tmp
            else:
                orders_heavy += orders_tmp
                torsions_heavy += torsions_tmp

        # sort in heavy-atom torsions followed by H-only torsions
        self.torsion_orders = orders_heavy + orders_H
        self.torsions = torsions_heavy + torsions_H

        if verbose:
            print('Torsions: ', self.torsions)

    def is_biaryl_ring(self,ring1,ring2):
        # first check if shares any atom
        for atm1 in ring1:
            for atm2 in ring2:
                if atm1 == atm2:
                    return False
        # Then figure out if any atm pair connected
        for atm1 in ring1:
            for atm2 in ring2:
                for bond in self.bonds:
                    if (atm1,atm2) != (bond.atm1,bond.atm2) and \
                            (atm1,atm2) != (bond.atm2,bond.atm1):
                        continue

                    #check if atm1,atm2 are connected by ring
                    is_connected_by_ring = False
                    for ring in self.rings:
                        if atm1 in ring and atm2 in ring:
                            is_connected_by_ring = True
                            break
                    if not is_connected_by_ring:
                        return (atm1,atm2)
        return False

    # add extra biaryl relation for ring-nonring (by going through hard-coded cases)
    def search_special_biaryl_ring(self,ring):
        biaryl_pivot_extra = []
        # [(A,B),..] : ring=A becomes special biaryl if any of ring=A-B connection is detected
        special_biaryl_to_ring = [(ACLASS_ID.index('CDp'),ACLASS_ID.index('OG2')),
                                  (ACLASS_ID.index('CDp'),ACLASS_ID.index('Oal')),
                                  (ACLASS_ID.index('CDp'),ACLASS_ID.index('Oad')),

                                  # ring-N=C
                                  (ACLASS_ID.index('Nin'),ACLASS_ID.index('CDp')),
                                  (ACLASS_ID.index('NG2'),ACLASS_ID.index('CDp')),
                                  (ACLASS_ID.index('Nam2'),ACLASS_ID.index('CDp')),
                                  (ACLASS_ID.index('Nin'),ACLASS_ID.index('CRp')),
                                  (ACLASS_ID.index('NG2'),ACLASS_ID.index('CRp')),
                                  (ACLASS_ID.index('Nam2'),ACLASS_ID.index('CRp')),

                                  # ring-N-SOn
                                  (ACLASS_ID.index('Nin'),ACLASS_ID.index('SG5')),
                                  (ACLASS_ID.index('NG2'),ACLASS_ID.index('SG5')),

                                  # ring-C=N*
                                  (ACLASS_ID.index('CDp'),ACLASS_ID.index('Nad')),
                                  (ACLASS_ID.index('CDp'),ACLASS_ID.index('Nam2')),
                                  (ACLASS_ID.index('CDp'),ACLASS_ID.index('NG2')),
                                  (ACLASS_ID.index('CDp'),ACLASS_ID.index('Nin')),

                                  (ACLASS_ID.index('CD1'),ACLASS_ID.index('CD1')),
                                  ]
        

        for a1 in ring:
            a2_is_special_biaryl = False
            for a2,order in self.atms[a1].bonds:
                a2class = self.atms[a2].aclass
                if a2 in self.atms_aro or a2 in self.atms_puckering: continue

                further_connected = False
                for a3,order in self.atms[a2].bonds:
                    #print self.atms[a1].name, self.atms[a2].name, self.atms[a3].name
                    if a3 == a1: continue
                    a3class = self.atms[a3].aclass
                    #print ring,self.atms[a1].name,self.atms[a2].name,self.atms[a3].name,
                    #ACLASS_ID[a2class], ACLASS_ID[a3class], (a2class,a3class) in special_biaryl_to_ring, len(self.atms[a3].bonds)
                    if (a2class,a3class) in special_biaryl_to_ring:
                        a2_is_special_biaryl = True
                        #break
                    if len(self.atms[a3].bonds) > 1:
                        further_connected = True
                #print self.atms[a1].name, self.atms[a2].name, a2_is_special_biaryl, further_connected
                if a2_is_special_biaryl and further_connected:
                    biaryl_pivot_extra.append((a1,a2))
                    break
        return biaryl_pivot_extra

    # this should be called when atomtypes are assigned
    def assign_rotable_torsions(self,verbose=False):
        self.rotable_torsion_id = []
        self.hapol_torsion_id = []
        self.hpol_torsion_type = [0 for k in self.torsions]

        for i,tor1 in enumerate(self.torsions[:-1]):
            for j,tor2 in enumerate(self.torsions[i+1:]):
                atms_in_list = []
                # use a funky way of checking members than recursive way...
                if (tor1[0] == tor2[0] and tor1[3] == tor2[3]) or (tor1[0] == tor2[3] and tor1[3] == tor2[0]): #6-member
                    atms_in_list = list(tor1)+[tor2[1],tor2[2]]
                elif (self.bond_order(tor1[0],tor1[3]) > 0): #4-member
                    atms_in_list = list(tor1)
                elif (tor1[0] == tor2[0] and tor1[3] == tor2[2]): #5-member
                    atms_in_list = list(tor1)+[tor2[1]]
                elif (tor1[0] == tor2[3] and tor1[3] == tor2[1]): #5-member
                    atms_in_list = list(tor1)+[tor2[2]]
                elif (tor1[0] == tor2[0] and self.bond_order(tor1[3],tor2[3]) > 0 ) or \
                     (tor1[3] == tor2[0] and self.bond_order(tor1[0],tor2[3]) > 0 ): #7-member
                    atm_in_list = list(tor1)+list(tor2[1:])
                elif (tor1[0] == tor2[3] and self.bond_order(tor1[3],tor2[0]) > 0 ) or \
                     (tor1[3] == tor2[3] and self.bond_order(tor1[0],tor2[0]) > 0 ):
                    atm_in_list = list(tor1)+list(tor2[:3])
                # skip >= 8-membered ring.. those are two flexible
                #elif (self.bond_order(tor1[0],tor2[0]) > 0 and self.bond_order(tor1[3],tor2[3]) > 0) or \
                #    (self.bond_order(tor1[0],tor2[3]) > 0 and self.bond_order(tor1[3],tor2[0]) > 0): #8-member
                #    atm_in_list = list(tor1)+list(tor2)
                else:
                    continue
                if atms_in_list == []: continue

                # make sure you are not doing it on duplicated atom list
                is_ring = True
                for atm in atms_in_list:
                    if atms_in_list.count(atm) > 1: 
                        is_ring = False
                        break
                if not is_ring: continue

                # check aromatic ring
                is_aro = True
                for atm in atms_in_list:
                    if self.atms[atm].hyb == 3 and (ACLASS_ID[self.atms[atm].aclass] not in ACLASS_ARO_SP3):
                        is_aro = False
                        break

                if is_aro:
                    atms_in_list_sort = copy.deepcopy(atms_in_list)
                    atms_in_list_sort.sort()
                    if atms_in_list_sort not in self.rings:
                        self.rings_aro.append(atms_in_list_sort)
                        self.rings.append(atms_in_list_sort)

                    for atm in atms_in_list:
                        if atm not in self.atms_aro:
                            self.atms_aro.append(atm)
                else: #otherwise puckering
                    for atm in atms_in_list:
                        if atm not in self.atms_puckering:
                            self.atms_puckering.append(atm)

                    # puckering also addes into rings (for the biaryl assignment...)
                    atms_in_list_sort = copy.deepcopy(atms_in_list)
                    atms_in_list_sort.sort()
                    if atms_in_list_sort not in self.rings:
                        self.rings.append(atms_in_list_sort)

        ###### search biaryl torsions
        self.biaryl_rings = []
        self.biaryl_pivots = []
        if len(self.rings) > 1:
            for i,ring1 in enumerate(self.rings[:-1]):
                for j,ring2 in enumerate(self.rings[i+1:]):
                    biaryl_pivot = self.is_biaryl_ring(ring1,ring2)
                    if biaryl_pivot:
                        self.biaryl_rings.append([i,i+j+1])
                        self.biaryl_pivots.append(biaryl_pivot)

            ## additional: non-ring connected to ring
            biaryl_pivot_extra = []
            for i,ring1 in enumerate(self.rings_aro):
                biaryl_pivot_extra += self.search_special_biaryl_ring(ring1)
            #print 'Added extra biaryl-pivots: ', 
            #for a1,a2 in biaryl_pivot_extra:
            #    print " (%s,%s)"%(self.atms[a1].name,self.atms[a2].name),
            #print
            self.biaryl_pivots += biaryl_pivot_extra
                
        # reassign atype
        for a1,a2 in self.biaryl_pivots:
            if self.atms[a1].aclass in [ACLASS_ID.index('CR'),ACLASS_ID.index('CRp'),
                                        ACLASS_ID.index('CD1'),ACLASS_ID.index('CD2'),ACLASS_ID.index('CD'),ACLASS_ID.index('CDp')]: #carbon
                self.atms[a1].aclass = ACLASS_ID.index('CRb')
            elif self.atms[a1].atype == ATYPES.index('N'):
                self.atms[a1].aclass = ACLASS_ID.index('NGb')
                
            if self.atms[a2].aclass in [ACLASS_ID.index('CR'),ACLASS_ID.index('CRp'),
                                        ACLASS_ID.index('CD1'),ACLASS_ID.index('CD2'),ACLASS_ID.index('CDp'),ACLASS_ID.index('CDp')]: #carbon
                self.atms[a2].aclass = ACLASS_ID.index('CRb')
            elif self.atms[a2].atype == ATYPES.index('N'):
                self.atms[a2].aclass = ACLASS_ID.index('NGb')

        if verbose:
            print('Atom_Puckering: ', [self.atms[atm].name for atm in self.atms_puckering])
            print('Atom_Aro: ', [self.atms[atm].name for atm in self.atms_aro])
            print('Rings: ', [[self.atms[atm].name for atm in ring] for ring in self.rings])
            if len(self.biaryl_rings) > 0 :
                print('BiarylRings:', self.biaryl_rings, ' BiarylAxes: ') 
                for a1,a2 in self.biaryl_pivots:
                    print((self.atms[a1].name,self.atms[a2].name))
                print()

        ###### hydrogens
        for i,torsion in enumerate(self.torsions):
            aclasses = [ACLASS_ID[self.atms[atm].aclass] for atm in torsion]
            # Apolar hydrogen torsions
            if (aclasses[0] in ACLASS_HAPOL) or (aclasses[3] in ACLASS_HAPOL):
                self.hapol_torsion_id.append(i)

            # PolarH torsions
            elif (aclasses[0] in ACLASS_HPOL) and (aclasses[3] not in ACLASS_HPOL):
                stem = self.atms[torsion[1]]
                if stem.hyb == 2:
                    if len(stem.bonds) < 3:
                        self.hpol_torsion_type[i] = 2
                        #self.rotable_torsion_id.append(i)
                elif stem.hyb == 3: #self.torsion_orders[i] == 1:
                    self.hpol_torsion_type[i] = 3
                    #self.rotable_torsion_id.append(i)

            elif (aclasses[3] in ACLASS_HPOL) and (aclasses[0] not in ACLASS_HPOL):
                stem = self.atms[torsion[2]]
                if stem.hyb == 2:
                    if len(stem.bonds) < 3:
                        self.hpol_torsion_type[i] = 2
                        #self.rotable_torsion_id.append(i)
                elif stem.hyb == 3: #self.torsion_orders[i] == 1:
                    self.hpol_torsion_type[i] = 3
                    #self.rotable_torsion_id.append(i)

    def bond_order(self,atm1,atm2):
        for bond in self.bonds:
            if (bond.atm1 == atm1 and bond.atm2 == atm2):
                return bond.order
            if (bond.atm1 == atm2 and bond.atm2 == atm1):
                return bond.order
        return 0

    def icoord(self,tup,FT):
        [i0,i1,i2,i3] = tup
        form = 'ICOOR_INTERNAL  %-4s %11.6f %11.6f %11.6f %-4s %-4s %-4s\n'
        xyz0 = self.xyz[i0]
        xyz1 = self.xyz[i1]
        xyz2 = self.xyz[i2]
        xyz3 = self.xyz[i3]
        if FT == 0:
            
            icoordstr =  form%(self.atms[i0].name,  0.0,  0.0,  0.0, 
                               self.atms[i0].name,self.atms[i1].name,self.atms[i2].name) # 1st
            #tup = (i0,i0,i1,i2)
        elif FT == 1:
            d = distance(xyz0,xyz1)
            icoordstr = form%(self.atms[i1].name,  0.0,180.0,    d, 
                              self.atms[i0].name,self.atms[i1].name,self.atms[i2].name) # 2nd
            #tup = (i1,i0,i1,i2)
        elif FT == 2:
            d = distance(xyz0,xyz1)
            a = 180.0-angle(xyz0,xyz1,xyz2)
            icoordstr = form%(self.atms[i0].name,  0.0,    a,    d, 
                              self.atms[i1].name,self.atms[i2].name,self.atms[i0].name) #3rd
            #tup = (i0,i1,i2,i0)
        else:
            d = distance(xyz0,xyz1)
            a = 180.0-angle(xyz0,xyz1,xyz2)
            t = dihedral(xyz0,xyz1,xyz2,xyz3)
            icoordstr = form%(self.atms[i0].name,    t,    a,    d, 
                              self.atms[i1].name,self.atms[i2].name,self.atms[i3].name) #3rd
            #tup = (i0,i1,i2,i3)
        return icoordstr

    def parse_crystinfo(self,l):
        self.crystinfo = ''

        words = l[:-1].split()
        cell = [float(word) for word in words[:6]]
        spacegrp = words[6]
        setting = int(words[7]) #not used

        SG_by_number = {
            #achiral
            '14':'P 1 21/c 1',
            '2':'P-1',
            '15':'C 1 2/c 1',
            '61':'Pbca',
            '33':'Pna21',
            '9':'C 1c 1',
            '60':'Pbcn',
            '29':'Pca21',
            '56':'Pccn',
            '231':'231',
            #chiral
            '19':'P 21 21 21',
            '4':'P 1 21 1',
            '1':'P 1',
            '5':'C 1 21',
            '18':'P 21 21 2',
        }

        # support when provided by number only for now
        if spacegrp in SG_by_number:
            self.crystinfo = 'CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f %s'%tuple(cell+[SG_by_number[spacegrp]])

    def report_pdbfile(self,outfile):
        out = open(outfile,'w')
        if self.crystinfo != '':
            out.write(self.crystinfo+'\n')

        form = 'HETATM%5d  %-4sLG1 X %3d    %8.3f%8.3f%8.3f  1.00 0.00\n'
        for i,atm in enumerate(self.atms):
            out.write(form%(i+1,atm.name,1,self.xyz[i][0],self.xyz[i][1],self.xyz[i][2]))
        out.close()

    def report_elec_cp_rep(self,outfile):
        out = file(outfile,'w')
        tip_index = [ACLASS_ID.index(a) for a in ['HS','HO','HN',
                                                  'Oal','Oad','OG2',
                                                  'Br','I','F','Cl',
                                                  'BrR','IR','FR','ClR']]
        for i,atm in enumerate(self.atms):
            if len(atm.bonds) > 1: continue

            if atm.aclass in tip_index:
                stem_atm = atm.bonds[0][0] 

                tip_with_same_stem = False
                for j,atm2 in enumerate(self.atms):
                    if i == j or len(atm2.bonds) > 1: continue
                    if atm2.bonds[0][0] == stem_atm:
                        tip_with_same_stem = True
                        break
                    if tip_with_same_stem: break
                if not tip_with_same_stem:
                    out.write('%s %4s %4s\n'%('LG1',#self.name,
                                              self.atms[stem_atm].name,atm.name)) #LG1 stem tip
        out.close()

    def define_icoord2(self,
                      Hapol_chi=False,
                      nbonded_chi=False,puckering_chi=True
                     ):
        FT_triple = []
        for i,atms in enumerate(self.angles):
            if self.nheavyatm >= 3 and \
               (self.atms[atms[0]].is_H or self.atms[atms[2]].is_H): continue
            if len(self.atms[atms[0]].bonds) == 1 and len(self.atms[atms[2]].bonds) == 1: continue
            if atms[0] == self.nbratom:
                FT_triple = list(atms)
                break
            elif atms[2] == self.nbratom:
                FT_triple = [atms[2],atms[1],atms[0]]
                break

        if FT_triple == []:
            for i,atms in enumerate(self.angles):
                if atms[1] == self.nbratom:
                    FT_triple = [atms[1],atms[0],atms[2]]

        # sorted for chi definition
        FTbase = {}
        FTbase[FT_triple[1]] = FT_triple[0]
        FTbase[FT_triple[2]] = FT_triple[1]

        # define atom roots spanning from FT_triple
        self.atms[FT_triple[2]].root = FT_triple[1]
        self.atms[FT_triple[1]].root = FT_triple[0]
        root_defined = copy.deepcopy(FT_triple)
        pivots_defined = []
        self.chiatms = []

        #WATCH OUT
        # look for chi on FT-defining atoms
        #FTpivot = (FT_triple[0],FT_triple[1])
        #for atms in self.torsions:
        #    if (atms[1],atms[2]) == FTpivot:
        #        #self.chiatms.append(list(atms))
                #pivots_defined.append(FTpivot)
        #        break
        #    elif (atms[2],atms[1]) == FTpivot:
                #self.chiatms.append([atms[3],atms[2],atms[1],atms[0]])
                #pivots_defined.append(FTpivot)
        #        break

        # define icoords, sort based on connectivity
        icoord_defs = []
        icoord_defined = copy.deepcopy(FT_triple) # is this different from root_defined?
        icoord_defs.append([FT_triple[0],FT_triple[1],FT_triple[2],FT_triple[0]])
        icoord_defs.append([FT_triple[0],FT_triple[1],FT_triple[2],FT_triple[0]])
        icoord_defs.append([FT_triple[2],FT_triple[1],FT_triple[0],FT_triple[2]])

        # first define for heavy atms
        while True:
            for bond in self.bonds:
                if self.atms[bond.atm1].is_H or self.atms[bond.atm2].is_H: continue
             
                chi_quad = []
                if bond.atm1 in root_defined and bond.atm2 not in root_defined:
                    root_defined.append(bond.atm2)
                    r0 = bond.atm2
                    self.atms[bond.atm2].root = bond.atm1
                    r1 = bond.atm1
                    r2 = self.atms[bond.atm1].root
                    #if broot > -1:
                    #    chi_quad = [self.atms[broot].root,broot,bond.atm1,bond.atm2]
                elif bond.atm2 in root_defined and bond.atm1 not in root_defined:
                    r0 = bond.atm1
                    root_defined.append(bond.atm1)
                    self.atms[bond.atm1].root = bond.atm2
                    r1 = bond.atm2
                    r2 = self.atms[bond.atm2].root
                    #if broot > -1: 
                    #    chi_quad = [self.atms[broot].root,broot,bond.atm2,bond.atm1]
                else:
                    continue

                if r0 in icoord_defined: continue
                # icoord definition
                if r1 == FT_triple[0]: #reverse root-direction
                    r2 = FT_triple[1]
                    r3 = FT_triple[2]
                elif r1 == FT_triple[1]: #attached to FT_triple[1]
                    #r2 = FT_triple[0]
                    r3 = FT_triple[2]
                elif r2 == FT_triple[0]:
                    r3 = FT_triple[1]
                else:
                    r3 = self.atms[r2].root
                icoord_defs.append([r0,r1,r2,r3])
                icoord_defined.append(r0)
                
                if len(root_defined) >= self.nheavyatm:
                    break
            if len(root_defined) >= self.nheavyatm:
                break

        # Then define roots for hydrogens
        for i,atm in enumerate(self.atms):
            if not atm.is_H: continue
            r1 = atm.bonds[0][0]
            r0 = i
            self.atms[i].root = r1
            r2 = self.atms[r1].root

            if r0 in icoord_defined: continue
            # icoord definition
            if r1 == FT_triple[0]: #reverse root-direction
                r2 = FT_triple[1]
                r3 = FT_triple[2]
            elif r1 == FT_triple[1]:
                r3 = FT_triple[2]
            elif r2 == FT_triple[0]:
                r3 = FT_triple[1]
            else:
                r3 = self.atms[r2].root
            icoord_defs.append([r0,r1,r2,r3])
            icoord_defined.append(r0)
                
        for icrd in icoord_defs:
            if -1 in icrd:
                print('WARNING: ICOORD not defined correctly for %s, %s!'%(self.name,self.atms[icrd[0]].name), icrd)
        #    print icrd

        # CHI
        self.get_chis1(FTbase,
                       Hapol_chi=Hapol_chi,
                       nbonded_chi=nbonded_chi,
                       puckering_chi=puckering_chi)

        return icoord_defs
        
    def get_chis2(self,
                  icoords,
                  Hapol_chi=False,
                  nbonded_chi=False,puckering_chi=True):

        #define chis
        num_H_confs = 1
        for i,atms in enumerate(self.chiatms):
            if self.hpol_torsion_type[i] == 2: num_H_confs *= 6
            elif self.hpol_torsion_type[i] == 3: num_H_confs *= 9
        if num_H_confs <= self.max_confs: self.chiextra = "1 20"
        else: self.chiextra = "0"
            
        self.chitypes = []
        for atms in copy.copy(self.chiatms):
            tipatm = self.atms[atms[3]]
            is_HPOLCHI = (ACLASS_ID[tipatm.aclass] in ACLASS_HPOL)
            is_HAPOLCHI = (ACLASS_ID[tipatm.aclass] in ACLASS_HAPOL)
            BIARYL_PIVOT = ((atms[1],atms[2]) in self.biaryl_pivots) or ((atms[2],atms[1]) in self.biaryl_pivots)
            if (not Hapol_chi) and is_HAPOLCHI:
                self.chiatms.remove(atms)
                continue
            if (not nbonded_chi) and (self.bond_order(atms[1],atms[2]) > 1) and not BIARYL_PIVOT:
                self.chiatms.remove(atms)
                continue 
            if (not puckering_chi) and (atms[1] in self.atms_puckering) and (atms[2] in self.atms_puckering):
                self.chiatms.remove(atms)
                continue

            root_hyb = self.atms[tipatm.root].hyb
            if is_HAPOLCHI: #sp3, aliphatic hydrogen
                self.chitypes.append('sp3H')
            elif is_HPOLCHI:
                if root_hyb == 3: #sp3
                    self.chitypes.append('sp3')
                #elif root_hyb  == 3: #sp2
                else:
                    self.chitypes.append('sp2')
                    if root_hyb not in [2,8,9]: 
                        print("Warning: hydrogen %s type unclear but assigned as sp2"%tipatm.name)
            elif (atms[1] in self.atms_puckering) and (atms[2] in self.atms_puckering):
                self.chitypes.append('pucker')
            else:
                self.chitypes.append('')

            
    def define_icoord1(self,
                      Hapol_chi=False,
                      nbonded_chi=False,puckering_chi=True
    ):
        ## ICOOR_INTERNAL
        #  let's make it as a crappy way for now...
        # define nbr atom first
        FT_triple = []
        for i,atms in enumerate(self.angles):
            if len(self.atms[atms[0]].bonds) == 1 and len(self.atms[atms[2]].bonds) == 1: continue
            if atms[0] == self.nbratom:
                FT_triple = list(atms)
                break
            elif atms[2] == self.nbratom:
                FT_triple = [atms[2],atms[1],atms[0]]
                break

        if FT_triple == []:
            for i,atms in enumerate(self.angles):
                if atms[1] == self.nbratom:
                    FT_triple = [atms[1],atms[0],atms[2]]
        
        if FT_triple == []:
            print('Unable to make FoldTree... exit')
            sys.exit()
        else:
            pass

        # sorted for chi definition
        FTbase = {}
        FTbase[FT_triple[1]] = FT_triple[0]
        FTbase[FT_triple[2]] = FT_triple[1]

        atoms_written = copy.deepcopy(FT_triple)

        nfail = 0
        chiatms = []
        icoord_defs = []
        icoord_defs.append([FT_triple[0],FT_triple[1],FT_triple[2],FT_triple[0]])
        icoord_defs.append([FT_triple[0],FT_triple[1],FT_triple[2],FT_triple[0]])
        icoord_defs.append([FT_triple[2],FT_triple[1],FT_triple[0],FT_triple[2]])

        while True:
            if len(self.torsions) == 0: break

            added=False
            for i,atms in enumerate(self.torsions):
                if (atms[3] not in atoms_written) and \
                   (atms[0] in atoms_written) and (atms[1] in atoms_written) and (atms[2] in atoms_written):
                    icoord_defs.append([atms[3],atms[2],atms[1],atms[0]])
                    chiatms.append(atms)
                    atoms_written.append(atms[3])
                    added=True
                    FTbase[atms[3]] = atms[2]
                    #print "Adding by torsion: ", self.atms[atms[3]].name
                    break
                
                if (atms[0] not in atoms_written) and \
                   (atms[1] in atoms_written) and (atms[2] in atoms_written) and (atms[3] in atoms_written):
                    icoord_defs.append([atms[0],atms[1],atms[2],atms[3]])
                    atoms_written.append(atms[0])
                    added=True
                    FTbase[atms[0]] = atms[1]
                    #print "Adding by torsion: ", self.atms[atms[0]].name
                    break

            # This shouldn't happen
            if not added:
                # add anything closest...
                for atm in self.atms_sorted_by_dist2nbratm:
                    if atm not in atoms_written:
                        print("Couldn't find torsion relation; adding by non-torsional relation: ", self.atms[atm].name)
                        icoord_defs.append([FT_triple[0],FT_triple[1],FT_triple[2],atm])
                        atoms_written.append(atm)
                        added=True
                        FTbase[atm] = FT_triple[2]

            if len(atoms_written) == len(self.atms):
                break
            if not added: nfail += 1
            
            if nfail > 10:
                print("Failure on assigning DOF!")
                sys.exit()

        # CHI
        self.get_chis1(FTbase,
                       Hapol_chi=Hapol_chi,
                       nbonded_chi=nbonded_chi,
                       puckering_chi=puckering_chi)

        return icoord_defs
        
    def get_chis1(self,
                  FTbase,
                  Hapol_chi=False,
                  nbonded_chi=False,puckering_chi=True
    ):

        self.chiatms = []
        self.chitypes = []
        self.chiextra = ""

        num_H_confs = 1
        covered = []
        for i,atms in enumerate(self.torsions):
            if (atms[1],atms[2]) in covered or (atms[2],atms[1]) in covered: continue
            if self.hpol_torsion_type[i] == 2: num_H_confs *= 6
            elif self.hpol_torsion_type[i] == 3: num_H_confs *= 9
            covered.append((atms[1],atms[2]))

        #print 'Total num_H_conf: ', num_H_confs
        if num_H_confs <= self.max_confs: self.chiextra = "1 20"
        else: self.chiextra = "0"
        
        covered = []

        # self.torsions are ordered such that heavyatom-only comes first followed by H-containing ones
        for i,atms in enumerate(self.torsions):
            # general rule
            if (atms[1],atms[2]) in covered or (atms[2],atms[1]) in covered: continue #avoid 
            if ((atms[1] in self.atms_aro) and (atms[2] in self.atms_aro)): continue #avoid part of aro
            # below are optional
            if (not Hapol_chi) and (i in self.hapol_torsion_id): continue #Apolar hydrogen chis
            if (not nbonded_chi) and (self.bond_order(atms[1],atms[2]) > 1): continue 
            if (not puckering_chi) and (atms[1] in self.atms_puckering) and (atms[2] in self.atms_puckering): continue
            #CHI

            #add chi only if connected through foldtree
            if not self.fold_tree_connection(atms):
                continue

            atm1 = self.atms[atms[1]]
            atm2 = self.atms[atms[2]]
            if atm1.root == atms[2]:
                atms_ordered = [atms[3],atms[2],atms[1],atms[0]] #last is tipatm
            else: # atm2.root = atms[1]
                atms_ordered = [atms[0],atms[1],atms[2],atms[3]] #last is tipatm

            # make sure sorted by index... this is Rosetta-specific thing
            #atms_ordered = copy.copy(atms)
            #if atms[0] in FTbase and atms[1] in FTbase and \
            #        FTbase[atms[0]] == atms[1] and FTbase[atms[1]] == atms[2]:
            #    atms_ordered = [atms[3],atms[2],atms[1],atms[0]]
            #elif atms[1] in FTbase and atms[2] in FTbase \
            #        and FTbase[atms[1]] == atms[0] and FTbase[atms[2]] == atms[1]:
            #    atms_ordered = [atms[0],atms[1],atms[2],atms[3]]

            self.chiatms.append(atms_ordered)

            if self.hpol_torsion_type[i] == 2: #sp2
                self.chitypes.append('sp2')
            elif self.hpol_torsion_type[i] == 3: #sp3
                self.chitypes.append('sp3')
            elif i in self.hapol_torsion_id: #sp3
                self.chitypes.append('sp3H')
            elif (atms[1] in self.atms_puckering) and (atms[2] in self.atms_puckering):
                self.chitypes.append('pucker')
            else:
                self.chitypes.append('')
            covered.append((atms[1],atms[2]))

    def setup_nbratom(self):
        com = [0.0 for k in range(3)]
        for crd in self.xyz[:self.nheavyatm]:
            for k in range(3):
                com[k] += crd[k]
        com = [x/self.nheavyatm for x in com]

        sortable = []
        for i,crd in enumerate(self.xyz[:self.nheavyatm+1]):
            sortable.append((distance(crd,com),i))
        sortable.sort()

        if len(self.atms) <= 2:
            self.nbratom = sortable[0][1]
        else:
            i = sortable[0][1]
            if len(self.atms[i].bonds) == 1:
                i = self.atms[i].bonds[0][0]
                if len(self.atms[i].bonds) == 1: # very unlikely
                    for d,i in sortable:
                        if len(self.atms[i].bonds) >= 2:
                            break
            self.nbratom = i
        #self.nbrradius = sortable[-1][0]+3.0
        self.nbrradius = (sortable[-1][0]+1.5)*2 #safe 
        self.atms_sorted_by_dist2com = [comp[1] for comp in sortable]

    def report_paramsfile(self,outfile,
                          report_Hapol_chi=False,
                          report_nbonded_chi=False,report_puckering_chi=True,
                          report_as_atype=False,report_as_mmtype=True,
                          verbose=True):
        out = open(outfile,'w')

        #resname = self.name
        resname = 'LG1'
        out.write('NAME %s\n'%resname)
        out.write('IO_STRING %s Z\n'%resname)
        out.write('TYPE LIGAND\n')
        out.write('AA UNK\n')

        # ATOM/BOND_TYPE
        for atom in self.atms:
            stdtype = MAP2STANDARD[ACLASS_ID[atom.aclass]] ##ATOM_TYPE
            gentype = ACLASS_ID[atom.aclass] #MM_ATOM_TYPE
            # logic since Jan 2017: use gentype as atype, mmtype as "X"
            #mmtype = gentype
            #atype = stdtype
            mmtype = "X"
            atype = gentype
            if report_as_atype:
                atype = gentype
            out.write('ATOM  %-4s %-4s %-4s %6.3f\n'%(atom.name,  #residue_io.cc, line 701 "atom_name(line.substr(5,4))"
                                                      atype,
                                                      mmtype,
                                                      atom.charge))
        for bond in self.bonds:
            if bond.order in [8,9]:
                border = 2
            else:
                border = bond.order
            atm1 = self.atms[bond.atm1]
            atm2 = self.atms[bond.atm2]
            out.write('BOND_TYPE %-4s %-4s %1d\n'%(atm1.name,atm2.name,border))

        #NBR setup
        self.setup_nbratom()
        out.write('NBR_ATOM %s\n'%self.atms[self.nbratom].name)
        out.write('NBR_RADIUS %8.5f\n'%self.nbrradius)

        # Get icoord and chiatm defs
        icoords = self.define_icoord2(report_Hapol_chi,
                                      report_nbonded_chi,
                                      report_puckering_chi)

        #ICOOR_INTERNAL    C4     0.000000    0.000000    0.000000   C4    C2    N1 
        #ICOOR_INTERNAL    C2     0.000000  180.000000    1.349077   C4    C2    N1 
        #ICOOR_INTERNAL    N1    -0.000001   52.329407    1.416959   C2    C4    N1
        for i,icoord_def in enumerate(icoords):
            out.write(self.icoord(icoord_def,i))

        for ichi,chiatms in enumerate(self.chiatms):
            #print ichi, chiatms, self.chitypes[ichi]
            atmnames = [self.atms[iatm].name for iatm in chiatms]
            atmnos = '%3d'*4%tuple(chiatms)
            comment = ''
            if self.chitypes[ichi] == 'pucker':
                comment = '#puckering'
            out.write('CHI %3d %4s %4s %4s %4s %s\n'%tuple([ichi+1]+atmnames+[comment])) #CHI 1  C2   C4   C5   C6
            
            #PROTON CHI
            if self.chitypes[ichi] == 'sp2':
                out.write("PROTON_CHI %i SAMPLES 2 0 180 EXTRA %s\n" % (ichi+1, self.chiextra))
            elif self.chitypes[ichi] == 'sp3':
                out.write("PROTON_CHI %i SAMPLES 3 60 -60 180 EXTRA %s\n" % (ichi+1, self.chiextra))
            elif self.chitypes[ichi] == 'sp3H':
                out.write("PROTON_CHI %i SAMPLES 3 60 -60 180 EXTRA %s\n" % (ichi+1, self.chiextra))

        out.close()

    def fold_tree_connection(self,atms):
        for i in range(len(atms)-1):
            a1,a2 = atms[i],atms[i+1]
            if (self.atms[a1].root != a2) and (self.atms[a2].root != a1):
                return False
        return True
        
    def define_chargegrps(self):
        self.chargegrps = []
        # first sort by abs(charge)
        sortable = []
        for i,atm in enumerate(self.atms):
            sortable.append((abs(atm.charge),i))
        sortable.sort()
        sortable.reverse()

        grped = []
        # let's make it conservative...
        for absq,iatm in sortable:
            if absq < 1.0: break

            is_grped = False
            for a2,border in self.atms[iatm].bonds:
                if a2 in grped:
                    is_grped = True
                    break

            if not is_grped:
                grp = [iatm]
                for a2,border in self.atms[iatm].bonds:
                    if abs(self.atms[a2].charge) > 0.3:
                        grp.append(a2)
                self.chargegrps.append(grp)
                grped += grp
                
        for i,atm in enumerate(self.atms):
            if i not in grped:
                grped.append(i)
                self.chargegrps.append([i])
        
    def report_grpdeffile(self,outfile):
        if self.chargegrps == []:
            self.define_chargegrps()

        out = file(outfile,'w')
        out.write('RESI LG1\n')
        for grp in self.chargegrps:
            grpstr = ' '.join([self.atms[i].name for i in grp])
            out.write('GROUP 0 %s\n'%(grpstr))
        out.close()
