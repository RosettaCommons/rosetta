#!/usr/bin/env python3
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
'''
Definitions or Basic classes used in mol2genparams.py or Molecule type.

Author: Hahnbeom Park and Frank DiMaio 
'''

import sys,os,math
import json
import operator
import importlib
from Types import *
from BasicClasses import OptionClass, AtomClass, BondClass, FunctionalGroupClass
from utils import distance,angle,dihedral

NUMBA_INSTALLED = (importlib.util.find_spec('numba'))
if NUMBA_INSTALLED: import AM1

import SetupTopology

class MoleculeClass:
    def __init__(self,mol2file,option=None):
        self.name = ''
        #self.resname = 'LG1' # goes to option
        self.mol2file = mol2file

        #read from mol2
        self.atms = []
        self.bonds = []
        self.xyz = []
        self.nheavyatm = 0
        self.max_confs = 5000 #copied from default of molfile_to_params
        self.crystinfo = ''

        #read from SetupTopology
        self.tree = None
        self.vatms = []
        self.angles = []
        self.torsions = []
        self.nbratom = -1
        self.nbrradius = 0.0
        self.atms_aro = []
        self.atms_puckering = []
        self.atms_ring = [] #including aromatic ring
        self.ATorder = []
        self.rings_aro = [] #only aromatic
        self.rings_pucker = [] #puckering
        self.rings = [] #aromatic + puckering
        self.biaryl_pivots = []
        self.biaryl_pivots_extra = []

        # chi-definition
        self.chiatms = []
        self.chitypes = []
        self.chiextra = ''

        # misc (optional)
        self.chargegrps = []
        self.funcgrps = []

        if option == None:
            self.option = OptionClass()
        else:
            self.option = option

        # self.bonds defined at read_mol2; angle,torsion at topology setup
        stat = self.read_mol2()
        if not stat:
            return

        SetupTopology.setup(self,option)
        # setup includes series of:
        #assign_bonds()
        #assign_hybridization_if_missing()
        #setup_nbratom()
        #assign_angles()
        #assign_torsions()
        #detect_rings()
        #define_icoord()
        #classifier.apply_to_molecule(mol)
        #classifier.assert_H(mol)
        #define_conjugation()
        #define_rotable_torsions()

        if self.option.opt.do_am1bcc_calc:
            self.calc_am1bcc_charge()

    def atom_index(self,atmname):
        for i,atm in enumerate(self.atms):
            if atm.name == atmname:
                return i
        return -1

    def read_mol2(self):
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
                if len(atms) > 0:
                    print("ERROR: More than one molecule defined in the mol2 file! Please use mol2 containing only one entry")
                    sys.exit()
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
                    order = 4 #amide as aromatic!
                elif words[3] == 'ar': #aromatic
                    order = 4 #or 2? 
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
        atomnames_read = {}
        for inew,iprv in enumerate(newindex):
            if self.option.verbose:
                print('Sorted atom index: %s %d -> %d'%(atms[iprv].name,iprv,inew))

            aname = atms[iprv].name
            if aname in atomnames_read:
                atomnames_read[aname].append(inew)
                newname = aname + '%d'%len(atomnames_read[aname])
                if self.option.verbose:
                    print('Renaming duplicate atom name %d %s -> %s'%(inew,aname,newname))
                atms[iprv].name = newname
            else:
                atomnames_read[aname] = [inew]
                
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

        return True

    def try_map_atype(self,atype):
        if atype.upper() == 'BR':
            return 'Br'
        elif atype.upper() == 'CL':
            return 'Cl'
        elif atype[0] in ATYPES_REG:
            return atype[0]
        else:
            for elem in SPECIAL_ATOMS:
                if atype == elem or atype.upper() == elem.upper():
                    return elem
            print('No such atom type: %s, while reading mol2file '%atype+self.mol2file+', line: ')
            return 'Null'

    def bond_order(self,atm1,atm2):
        for bond in self.bonds:
            if (bond.atm1 == atm1 and bond.atm2 == atm2):
                return bond.order
            if (bond.atm1 == atm2 and bond.atm2 == atm1):
                return bond.order
        return 0

    def bond_conjugated(self,atm1,atm2):
        for bond in self.bonds:
            if (bond.atm1 == atm1 and bond.atm2 == atm2):
                return bond.is_conjugated
            if (bond.atm1 == atm2 and bond.atm2 == atm1):
                return bond.is_conjugated
        return False

    # 4 for aro, 1 for puckering, 0 for none
    def in_same_ring(self,atm1,atm2):
        for ring in self.rings:
            if ring.has([atm1,atm2]):
                return ring.type
        return 0 # not in same ring

    def report_icoord(self,outstream):
        icoordcont = []
        form = "ICOOR_INTERNAL  %-4s %11.6f %11.6f %11.6f %-4s %-4s %-4s\n"
        for i,iatm in enumerate(self.ATorder):
            atm = self.atms[iatm]
            len_i,ang_i,dih_i=0.0,180.0,0.0
            #print ("ROOT:", iatm, atm.root, atm.groot )
            if (i > 0):
                len_i = distance(self.xyz[iatm],self.xyz[atm.root])
                if (i > 1):
                    ang_i = angle(self.xyz[iatm],self.xyz[atm.root],self.xyz[atm.groot[0]])
                    if (i > 2):
                        dih_i = dihedral(self.xyz[iatm],self.xyz[atm.root],
                                         self.xyz[atm.groot[0]],self.xyz[atm.groot[1]])

            l = form%(atm.name,dih_i,180.0-ang_i,len_i,
                      self.atms[atm.root].name,self.atms[atm.groot[0]].name,self.atms[atm.groot[1]].name)
            outstream.write(l)

        # virtual atms
        if self.option.opt.report_puckering_chi:
            for i,atm in enumerate(self.vatms):
                ring = self.rings[atm.ring_index]
                if ring.type != 1: continue
                
                ivrt = atm.vrt_i
                #print ("VROOT:", ivrt, atm.root, atm.groot)
                len_i = distance(self.xyz[ivrt],self.xyz[atm.root])
                ang_i = angle(self.xyz[ivrt],self.xyz[atm.root],self.xyz[atm.groot[0]])
                dih_i = dihedral(self.xyz[ivrt],self.xyz[atm.root],
                                 self.xyz[atm.groot[0]],self.xyz[atm.groot[1]])

                l = form%(atm.name,dih_i,180.0-ang_i,len_i,
                          self.atms[atm.root].name,self.atms[atm.groot[0]].name,self.atms[atm.groot[1]].name)
                outstream.write(l)
                        
    def parse_crystinfo(self,l):
        self.crystinfo = ''

        words = l[:-1].split()
        cell = [float(word) for word in words[:6]]
        spacegrp = words[6]
        setting = int(words[7]) #not used

        direc = os.path.dirname( os.path.abspath(__file__) )
        spacegroupfn = os.path.join(direc, "spacegroup.json")
        if os.path.exists(spacegroupfn):            
            with open( spacegroupfn ) as f:
                SG_by_number = json.load(f)
        else:
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
        else:
            raise RuntimeError("No space group %s"%spacegrp)
            
    def report_pdbfile(self,outfile):
        out = open(outfile,'w')
        if self.crystinfo != '':
            out.write(self.crystinfo+'\n')

        form = 'HETATM%5d  %-4s%3s X %3d    %8.3f%8.3f%8.3f  1.00 0.00\n'
        for i,atm in enumerate(self.atms):
            out.write(form%(i+1,atm.name,self.option.get_resname(),1,
                            self.xyz[i][0],self.xyz[i][1],self.xyz[i][2]))
        out.close()

    def report_elec_cp_rep(self,outfile):
        out = open(outfile,'w')
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
                    out.write('%s %4s %4s\n'%(self.option.get_resname(),
                                              self.atms[stem_atm].name,atm.name)) #LG1 stem tip
        out.close()

    def report_functional_grps(self,out):
        for grp in self.funcgrps:
            grp.show(out)

    def report_paramsfile(self,outfile):
        out = open(outfile,'w')

        # start writing
        resname = self.option.get_resname()
        out.write('NAME %s\n'%resname)
        out.write('IO_STRING %s Z\n'%resname)
        out.write('TYPE LIGAND\n')
        out.write('AA UNK\n')

        if self.option.opt.report_puckering_chi and self.vatms != []:
            out.write('PROPERTIES CYCLIC\n')

        # ATOM
        for atom in self.atms:
            gentype = ACLASS_ID[atom.aclass] 
            out.write('ATOM %-4s %-4s %-4s %6.3f\n'%(atom.name,
                                                     gentype,'X',
                                                     atom.charge))

        # ATOMS-VIRTUAL
        if self.option.opt.report_puckering_chi:
            for atom in self.vatms:
                ring = self.rings[atom.ring_index]
                if ring.type != 1: continue 
                out.write('ATOM %-4s %-4s %-4s %6.3f\n'%(atom.name,
                                                         'VIRT','VIRT',
                                                         atom.charge))
            
        #BOND_TYPE
        for bond in self.bonds:
            if bond.order not in BOND_ORDERS:
                sys.exit('Unknown bond order of %d!'%bond.order)
            border = bond.order_in_params()
            atm1 = self.atms[bond.atm1]
            atm2 = self.atms[bond.atm2]
            #print (self.atms[bond.atm1].name,self.atms[bond.atm2].name,border,bond.is_conjugated)

            # Assign "RING" info for atom pair in a puckering ring (== not aromatic) 
            extra = ''
            if self.in_same_ring(bond.atm1,bond.atm2) > 0:
                extra = ' RING'
            out.write('BOND_TYPE %-4s %-4s %1d%s\n'%(atm1.name,atm2.name,border,extra))

        #VIRTUAL BONDS
        if self.option.opt.report_puckering_chi:
            for rnum,ring in enumerate(self.rings):
                if ring.type != 1: continue

                (k,l) = ring.cut_bond
                ktag = "V%du"%rnum
                ltag = "V%dl"%rnum
                out.write('BOND_TYPE %-4s %-4s %1d\n'%(self.atms[k].name,ltag,1))
                out.write('BOND_TYPE %-4s %-4s %1d\n'%(self.atms[l].name,ktag,1))
                out.write("VIRTUAL_SHADOW %-4s %-4s\n"%(self.atms[k].name,ktag))
                out.write("VIRTUAL_SHADOW %-4s %-4s\n"%(self.atms[l].name,ltag))

        # CUT_BOND
        cut_bonds_reported = []
        for rnum,ring in enumerate(self.rings):
            (i,j) = ring.cut_bond
            if (i,j) in cut_bonds_reported or (j,i) in cut_bonds_reported:
                continue
            cut_bonds_reported.append((i,j))
            out.write('CUT_BOND %4s %4s\n'%(self.atms[i].name,self.atms[j].name))

        # NBR
        out.write('NBR_ATOM %s\n'%self.atms[self.nbratom].name)
        out.write('NBR_RADIUS %8.5f\n'%self.nbrradius)

        # Report ICOORD following bond definitions
        #C4: center of angle
        #ICOOR_INTERNAL    C4     0.000000    0.000000    0.000000   C4    C2    N1 
        #ICOOR_INTERNAL    C2     0.000000  180.000000    1.349077   C4    C2    N1 
        #ICOOR_INTERNAL    N1    -0.000001   52.329407    1.416959   C4    C2    N1
        self.report_icoord(out)

        # RING
        ringchis_counted = []
        if self.option.opt.report_puckering_chi:
            nuctr = 0
            ringctr = 0
            for rnum,ring in enumerate(self.rings):
                ringchis_counted.append(ring)
                
                if ring.type != 1:
                    if ring.type in [2,3]:
                        out.write('#RING %2d %3d'%(ring.type,ring.natms)+\
                                  ' %4s'*ring.natms%tuple([self.atms[i].name for i in ring.atms])+'\n')
                    continue

                ringctr += 1
                ktag = "V%du"%(rnum)
                ltag = "V%dl"%(rnum)

                l = "ADD_RING %d"%(ringctr)
                for i in ring.atms:
                    l += " "+self.atms[i].name
                out.write(l+"\n")

                for i in range(ring.natms-1):
                    tag1 = self.atms[ring.atms[i-1]].name
                    if (i-1<0):
                        tag1 = ktag
                    tag2 = self.atms[ring.atms[i]].name
                    tag3 = self.atms[ring.atms[i+1]].name
                    tag4 = self.atms[ring.atms[(i+2)%ring.natms]].name
                    if (i+2>=ring.natms):
                        tag4 = ltag
                    nuctr+=1
                    l = "NU %2d %-4s %-4s %-4s %-4s\n"%(nuctr,tag1,tag2,tag3,tag4)
                    out.write(l)

        # CHI
        ichi = 0
        for i,chiatms in enumerate(self.chiatms):
            atmnames = [self.atms[iatm].name for iatm in chiatms]
            atmnos = '%3d'*4%tuple(chiatms)

            # skip if already defined as ringchis
            comment = ''
            counted = False
            for ring in ringchis_counted:
                if ring.has(chiatms[1:3]):
                    counted = True
                    break
            if counted:
                continue
            #else:
            #    comment = '#puckering'
            
            if (chiatms[1],chiatms[2]) in self.biaryl_pivots or \
               (chiatms[2],chiatms[1]) in self.biaryl_pivots:
               comment = "#biaryl"
                
            out.write('CHI %3d %4s %4s %4s %4s %s\n'%tuple([ichi+1]+atmnames+[comment])) #CHI 1  C2   C4   C5   C6
            
            #PROTON CHI
            if self.chitypes[ichi] == 'sp2':
                out.write("PROTON_CHI %i SAMPLES 2 0 180 EXTRA %s\n" % (ichi+1, self.chiextra))
            elif self.chitypes[ichi] == 'sp3':
                out.write("PROTON_CHI %i SAMPLES 3 60 -60 180 EXTRA %s\n" % (ichi+1, self.chiextra))
            elif self.chitypes[ichi] == 'sp3H':
                out.write("PROTON_CHI %i SAMPLES 3 60 -60 180 EXTRA %s\n" % (ichi+1, self.chiextra))
            # unused long-ring: commented out from params
            ichi += 1

        out.close()

    #should I merge define_chis1 & assign_rotable_torsions?
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

        out = open(outfile,'w')
        out.write('RESI LG1\n')
        for grp in self.chargegrps:
            grpstr = ' '.join([self.atms[i].name for i in grp])
            out.write('GROUP 0 %s\n'%(grpstr))
        out.close()

    def calc_am1bcc_charge(self):
        if not NUMBA_INSTALLED:
            print( "AM1BCC calculation requested but numba library not installed! skip." )
            return
        a = AM1.AM1(open(self.mol2file),streamtype="mol2")
        a.runLBFGS()
        Z = a.get_charges()
        
        for i,atmname in enumerate(a.mol.atms):
            iatm = self.atom_index(atmname)
            if iatm < 0:
                print( "Warning while converting AM1 charge: atom index unknown! %s"%(atmname) )
            self.atms[iatm].charge = Z[i] #overwrite
            
