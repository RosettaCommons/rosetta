import sys,os
from utils import dihedral
from Types import *
MYFILE = os.path.abspath(__file__)
direc = MYFILE.replace(MYFILE.split('/')[-1],'')

# Wildcard for torsion assignment but not used now, replaced to gen_torsion code
WILDCARDS = {'C*':['CS','CS1','CS2','CS3','CD','CD1','CD2','CR','CT','CSp','CDp','CRp','CTp'],
             'H*':['HC','HR','HO','HN','HS','HG'],
             'N*':['Nam','Nam2','Nad','Nad3','Nin','Nim','Ngu1','Ngu2','NG1','NG2','NG21','NG22','NG3','NGb'],
             'O*':['Ohx','Oet','Oal','Oad','Oat','Ofu','OG2','OG3','Ont'],
             'P*':['PG3','PG5'],
             'S*':['Sth','Ssl','SG2','SG3','SG5'],
             'X':copy.copy(ACLASS_ID)
             }
WILDCARD_ATOMS = list(WILDCARDS.keys())
WILDCARD_ATOMS.sort()

WILDCARD_INDEX = {}
for i,atm in enumerate(ACLASS_ID):
    WILDCARD_INDEX[i] = []

for i,atm in enumerate(ACLASS_ID):
    for j,wc in enumerate(WILDCARD_ATOMS):
        if atm in WILDCARDS[wc]:
            WILDCARD_INDEX[i].append(NCLASS+j)

#for i,atm in enumerate(ACLASS_ID):
#    print 'WC setup: ', i, atm, WILDCARD_INDEX[i]

#### Torsion grouping setup -- not used
#unshared: Ohx, HX, Oad, Ofu, ...
SHARE_TORSIONGROUP = [['Null'], #keep this!!
                      ['CS','CS1','CS2','CS3','CSp'],
                      ['CD','CD1','CD2','CDp'],
                      ['CR','CRp','CRb'],
                      ['CT','CTp'],
                      ['OG2','Oal','Ofu'], # any?
                      ['OG3','Oet','OG31'],
                      ['NG2','NG21','NG22','Nin','Nim','Nad','Nad3','Ngu1','Ngu2','Nam2','NGb'],
                      ['NG3','Nam','Nam2'],
                      ['SG3','Ssl','SG5'],
                      ['PG3','PG5'],
                      #['HG','HC','HR','HO','HN','HS'],
                  ]

TORSION_GROUP = {}
ungrouped = 0
TORSION_GROUP_REPS = {}

SIMPLE_TORS = True
if SIMPLE_TORS:
    for iatm,atm in enumerate(ACLASS_ID):
        found_group = False
        for igrp,group in enumerate(SHARE_TORSIONGROUP):
            if atm in group:
                found_group = True
                break
        if found_group:
            TORSION_GROUP[iatm] = igrp
            TORSION_GROUP_REPS[igrp] = SHARE_TORSIONGROUP[igrp][0] #string
        else:
            igrp = len(SHARE_TORSIONGROUP) + ungrouped
            TORSION_GROUP[iatm] = igrp
            TORSION_GROUP_REPS[igrp] = atm #string
            ungrouped += 1
else:
    for iatm,atm in enumerate(ACLASS_ID+WILDCARD_ATOMS):
        TORSION_GROUP[iatm] = iatm
        TORSION_GROUP_REPS[iatm] = atm #string

#### Torsion Pivot setup
aclasses_in_torspivot = [['CS','CS1','CS2','CS3','CSp'], #0: CS
                         ['CD','CD1','CD2','CDp'], #1: CD
                         ['CR','CRp','CRb'], #2: CR
                         ['CS','CSp','CD','CDp','CR','CRp','CT','CTp'], #3:C*
                         ['NG3','Nam','Nam2'], #4: NG3
                         ['NG2','NG21','NG22','Nin','Nim','Nad','Nad3','Ngu1','Ngu2','Nam2','NGb'], #5: NG2
                         ['OG2','OG3','Oal','Oet','Ohx','Ofu','OG31'], #6: O*
                         ['SG2','SG3','SG5','Ssl'], #7: S*
                         ['PG3','PG5'], #8: P*
                         copy.copy(ACLASS_ID), #9: X
]

#1(CS) -> [0,3,9], ...
TORSION_PIVOT_ID = {}
for iatm in TORSION_GROUP:
    TORSION_PIVOT_ID[iatm] = []
    for ipivot,aclasses in enumerate(aclasses_in_torspivot):
        #if iatm not in ACLASS_ID: continue #wildcard
        if ACLASS_ID[iatm] in aclasses:
            TORSION_PIVOT_ID[iatm].append(ipivot)

#for iatm in TORSION_PIVOT_ID:
#    print ACLASS_ID[iatm], iatm, TORSION_PIVOT_ID[iatm]

class TorsionParam:
    def __init__(self,name,atmtypes,ks):
        self.name = name
        self.index = int('%02d'*4%atmtypes)
        self.ks = ks

class TorsionDatabase:
    def __init__(self):
        self.db = {}
        for k in range(100):
            self.db[k] = []

    def add(self,torsclass,name,atmtypes,ks):
        self.db[torsclass].append(TorsionParam(name,atmtypes,ks))

    def search(self,torsclass,aclasses):
        myindex = int('%02d'*4%aclasses)

        # first search for explicit type
        for i,param in enumerate(self.db[torsclass]):
            if param.index == myindex:
                return i

        # If fails, search for wild types
        wc1s = [aclasses[0]]+WILDCARD_INDEX[aclasses[0]]
        wc2s = [aclasses[3]]+WILDCARD_INDEX[aclasses[3]]

        pivots_avail = []
        for piv1 in [aclasses[1]]+WILDCARD_INDEX[aclasses[1]]:
            for piv2 in [aclasses[2]]+WILDCARD_INDEX[aclasses[2]]:
                pivots_avail.append((piv1,piv2))

        #print 'searching for %d*%d*%d = %d combinations'%(len(wc1s),len(wc2s),len(pivots_avail),
        #                                                  len(wc1s)*len(wc2s)*len(pivots_avail) )

        for i,param in enumerate(self.db[torsclass]):
            for pivot in pivots_avail:
                for wc1 in wc1s:
                    for wc2 in wc2s:
                        myindex = int('%02d'*4%(wc1,pivot[0],pivot[1],wc2))
                        #if torsclass == 56:
                        #    print i,param.index,myindex
                        if param.index == myindex:
                            return i
        # If all fails
        return -1

    def get_param(self,torsclass,torsid):
        return self.db[torsclass][torsid]

    def size(self):
        n = 0
        for torsclass in list(self.db.keys()):
            n += len(self.db[torsclass])
        return n

class TorsionAssigner:
    def __init__(self,dbfile=""):
        self.read_db(dbfile)
        
    def read_db(self,dbfile):
        self.torsdb = TorsionDatabase()
        TORSATOMS = ACLASS_ID+WILDCARD_ATOMS

        readtors = False
        readspecial = False
        
        for l in open(dbfile):
            if l.startswith('SPECIAL_TORSION'):
                readspecial = True
                continue
            if l.startswith('TORSION'):
                readtors = True
                readspecial = False
                continue
            if l.startswith('IMPROPER'): break
            if l.startswith('#'): continue

            if not (readtors or readspecial): continue

            words = l[:-1].split()
            if len(words) < 8: continue
            if readspecial:
                torsclass = 0
                a1,a2,a3,a4 = words[:4]
                k1,k2,k3 = [float(k) for k in words[5:8]]
            else:
                torsclass = int(words[0])
                a1,a2,a3,a4 = words[1:5]
                k1,k2,k3 = [float(k) for k in words[6:9]]

            indices = [TORSATOMS.index(a) for a in [a1,a2,a3,a4]]
            # make sure that all the indices are in torsion_group
            for i,iatm in enumerate(indices):
                if (iatm not in TORSION_GROUP) and (TORSATOMS[iatm] not in WILDCARD_ATOMS):
                    print('ERROR: Atom %s not defined as TORSION_GROUPS nor WILDCARD_ATOMS! '%(words[:4][i]))
                    sys.exit()

            # use torsion group instead
            if indices[1] > indices[2]:
                indices = [indices[3],indices[2],indices[1],indices[0]]
            elif indices[1] == indices[2]:
                if indices[0] > indices[3]:
                    indices = [indices[3],indices[2],indices[1],indices[0]]
            name = '-'.join([TORSATOMS[i] for i in indices])

            self.torsdb.add(torsclass,name,tuple(indices),(k1,k2,k3))
            #print 'ADD DB:', name, '%02d'%torsclass, self.torsdb.db[torsclass][-1].index
        
    def get_torsionclass(self,i1,i2):
        classes = []
        # make all possible including C*,X

        for j1 in TORSION_PIVOT_ID[i1]:
            for j2 in TORSION_PIVOT_ID[i2]:
                if j1 <= j2:
                    classes.append(j1*10+j2)
        return classes

    def assign(self,molecule,debug=False):
        molecule.torsionparams = []
        for itor,(a1,a2,a3,a4) in enumerate(molecule.torsions):
            atm1 = ACLASS_ID[molecule.atms[a1].aclass]
            atm2 = ACLASS_ID[molecule.atms[a2].aclass]
            atm3 = ACLASS_ID[molecule.atms[a3].aclass]
            atm4 = ACLASS_ID[molecule.atms[a4].aclass]
            #print 'Work %3s %3s %3s %3s'%(atm1, atm2, atm3, atm4)

            if (a2,a3) in molecule.biaryl_pivots or (a3,a2) in molecule.biaryl_pivots:
                if debug:
                    dihe = dihedral(molecule.xyz[a1],molecule.xyz[a2],
                                    molecule.xyz[a3],molecule.xyz[a4])

                    if 180.0 - abs(dihe) > 30.0 and abs(dihe) > 30.0:
                        print("Warning, biaryl torsion: ", itor, a1,a2,a3,a4, dihe)
                

            aclasses = [molecule.atms[a1].aclass,
                        molecule.atms[a2].aclass,
                        molecule.atms[a3].aclass,
                        molecule.atms[a4].aclass]
            
            # switch to torsion representatives
            # this is pretty messy... torsion_group_id -> group_rep_atom_name -> index of the rep_atom
            tmp = copy.deepcopy(aclasses)
            aclasses = [ ACLASS_ID.index(TORSION_GROUP_REPS[TORSION_GROUP[aclass]]) for aclass in aclasses]

            if aclasses[1] > aclasses[2] or \
               ((aclasses[1] == aclasses[2]) and (aclasses[3] < aclasses[0])):
                aclasses = [aclasses[3],aclasses[2],aclasses[1],aclasses[0]]

            torsclasses = self.get_torsionclass(aclasses[1],aclasses[2])
            #print ACLASS_ID[aclasses[1]], ACLASS_ID[aclasses[2]], torsclasses

            for i,torsclass in enumerate(torsclasses):
                torsid = self.torsdb.search(torsclass, tuple(aclasses))
                if torsid >= 0:
                    break

            if torsid < 0:
                anames = [ACLASS_ID[i] for i in aclasses]
                if debug:
                    print('Failed on: %s,'%molecule.name+' %s-%s-%s-%s'%tuple(anames)+', TorsionClasses: ', torsclasses)
                molecule.torsionparams.append((torsclass,-1))
            else:
                torsparam = self.torsdb.get_param(torsclass,torsid)
                if debug:
                    index = int('%02d'*4%tuple(aclasses))
                    atmnos = '%3d'*4%(a1,a2,a3,a4)
                    form = '%s Tors %3d, Assign %s (class %2d, id %3d %8d:%3s %3s %3s %3s) -> %-15s'
                    print(form%(molecule.name,itor,atmnos,torsclass,torsid,index,
                                atm1,atm2,atm3,atm4,torsparam.name))
                molecule.torsionparams.append((torsclass,torsid))
                if torsparam.name == 'X-X-X-X':
                    print(atm1, atm2, atm3, atm4)
