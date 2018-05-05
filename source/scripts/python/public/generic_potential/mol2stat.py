import sys,copy

from utils import dihedral, angle, dat2distribution, distance
from Types import *
from Molecule import MoleculeClass
from AtomTypeClassifier import AtomTypeClassifier
from TorsionAssigner import TorsionAssigner, TORSION_GROUP

class StatClass:
    def __init__(self):
        self.atype_counts = [0 for a in ACLASS_ID]
        self.torsions = []
        self.torsion_counts = []
        self.bonds = []
        self.bond_dstr = []
        self.angles = []
        self.angle_dstr = []
        self.torsion_params = []
        self.torsion_param_counts = []
        self.torsion_dstr = []

        self.MINTORS = -180.0
        self.MAXTORS = 180.0
        self.TORSBIN = 15.0
        self.NTORSBIN = int((self.MAXTORS-self.MINTORS)/self.TORSBIN)

        self.MINANG = 90.0
        self.MAXANG = 135.0
        self.ANGBIN = 2.5
        self.NANGBIN = int((self.MAXANG-self.MINANG)/self.ANGBIN)

        self.MINBOND = 0.8
        self.MAXBOND = 2.2
        self.BONDBIN = 0.05
        self.NBONDBIN = int((self.MAXBOND-self.MINBOND)/self.BONDBIN)

    def init_param(self,torsdb):
        self.torsion_params = []
        self.torsion_param_counts = []
        for torsclass in torsdb:
            for i,param in enumerate(torsdb[torsclass]):
                self.torsion_params.append((torsclass,i))
                self.torsion_param_counts.append(0)
                self.torsion_dstr.append([])
        
    def append_stat(self,molecule):
        atms = molecule.atms
        for atm in atms:
            self.atype_counts[atm.aclass] += 1

        for torsclass,itors in molecule.torsionparams:
            i = self.torsion_params.index((torsclass,itors))
            self.torsion_param_counts[i] += 1

        for i,bond in enumerate(molecule.bonds):
            a1 = bond.atm1
            a2 = bond.atm2
            xyz1 = molecule.xyz[a1]
            xyz2 = molecule.xyz[a2]
            d = distance(xyz1,xyz2)

            bondtype = [molecule.atms[a].aclass for a in [a1,a2]]
            if bondtype[0] > bondtype[1]:
                bondtype = [bondtype[1],bondtype[0]]
                
            if bondtype not in self.bonds:
                self.bonds.append(bondtype)
                self.bond_dstr.append([])

            ibond = self.bonds.index(bondtype)
            self.bond_dstr[ibond].append(d)

        for i,(a1,a2,a3) in enumerate(molecule.angles):
            xyz1 = molecule.xyz[a1]
            xyz2 = molecule.xyz[a2]
            xyz3 = molecule.xyz[a3]
            ang = angle(xyz1,xyz2,xyz3)

            angtype = [molecule.atms[a].aclass for a in [a1,a2,a3]]
            if angtype[0] > angtype[2]:
                angtype = [angtype[2],angtype[1],angtype[0]]
                
            if angtype not in self.angles:
                self.angles.append(angtype)
                self.angle_dstr.append([])

            iang = self.angles.index(angtype)
            self.angle_dstr[iang].append(ang)
            
        for i,(a1,a2,a3,a4) in enumerate(molecule.torsions):
            xyz1 = molecule.xyz[a1]
            xyz2 = molecule.xyz[a2]
            xyz3 = molecule.xyz[a3]
            xyz4 = molecule.xyz[a4]
            torsion = dihedral(xyz1,xyz2,xyz3,xyz4)

            torsclass,itors = molecule.torsionparams[i]
            itors_tot = self.torsion_params.index((torsclass,itors))
            self.torsion_dstr[itors_tot].append(torsion)

            torstype = (TORSION_GROUP[atms[a1].aclass],
                        TORSION_GROUP[atms[a2].aclass],
                        TORSION_GROUP[atms[a3].aclass],
                        TORSION_GROUP[atms[a4].aclass])

            if torstype[0] == 0 or torstype[1] == 0 or torstype[2] == 0 or torstype[3] == 0:
                continue 

            # make sure non-redundant
            # by making first one having lower index than the final
            if torstype[0] > torstype[3]:
                torstype_cp = copy.copy(torstype)
                torstype = (torstype_cp[3],
                            torstype_cp[2],
                            torstype_cp[1],
                            torstype_cp[0])
                
            if torstype not in self.torsions:
                self.torsions.append(torstype)
                self.torsion_counts.append(0)

            itors = self.torsions.index(torstype)
            self.torsion_counts[itors] += 1

    def report_atypes(self):
        for i,atype in enumerate(ACLASS_ID):
            print('%2d %-4s %5d'%(i,atype,self.atype_counts[i]))

    def report_bond_params(self):
        for i,bond in enumerate(self.bonds):
            bondname = '.'.join([ACLASS_ID[val] for val in bond])
            if 'Null' in bondname: continue
            l='%3d %-15s %5d %6.3f %6.3f|'%(i,bondname,
                                            len(self.bond_dstr[i]),
                                            0.0, 0.0)

            if len(self.bond_dstr[i]) > 0:
                try:
                    dstr = dat2distribution(self.bond_dstr[i],
                                            self.MINBOND,self.MAXBOND,self.BONDBIN)
                    dstr = [val*100.0 for val in dstr]
                    l += ' %4.1f'*self.NBONDBIN%tuple(dstr)
                except:
                    pass
            print(l)

    def report_angle_params(self):
        for i,ang in enumerate(self.angles):
            angname = '.'.join([ACLASS_ID[val] for val in ang])
            if 'Null' in angname: continue
            l='%3d %-15s %5d %6.3f %6.3f %6.3f|'%(i,angname,
                                                  len(self.angle_dstr[i]),
                                                  0.0, 0.0, 0.0)

            if len(self.angle_dstr[i]) > 0:
                try:
                    dstr = dat2distribution(self.angle_dstr[i],
                                            self.MINANG,self.MAXANG,self.ANGBIN)
                    dstr = [val*100.0 for val in dstr]
                    l += ' %4.1f'*self.NANGBIN%tuple(dstr)
                except:
                    pass
            print(l)

    def report_torsion_params(self,torsdb):
        for i,(torsclass,itors) in enumerate(self.torsion_params):
            p = torsdb[torsclass][itors]
            l='%2d %3d %-15s %5d %6.3f %6.3f %6.3f|'%(torsclass,itors,p.name,
                                                      self.torsion_param_counts[i],
                                                      p.ks[0], p.ks[1], p.ks[2])

            if len(self.torsion_dstr[i]) > 0:
                try:
                    dstr = dat2distribution(self.torsion_dstr[i],
                                            self.MINTORS,self.MAXTORS,self.TORSBIN)
                    dstr = [val*100.0 for val in dstr]
                    l += ' %4.1f'*self.NTORSBIN%tuple(dstr)
                except:
                    pass
            print(l)
        
    def report_torsions(self):
        # sort by abundance before reporting
        sortable = [[self.torsion_counts[i],torsion] for i,torsion in enumerate(self.torsions)]
        sortable.sort()
        sortable.reverse()
        self.torsions = [comp[1] for comp in sortable]
        self.torsion_counts = [comp[0] for comp in sortable]

        form = '%2d %4s %4s %4s %4s %5d'
        for i,(t1,t2,t3,t4) in enumerate(self.torsions):
            print(form%(i,
                        TORSION_GROUP_REPS[t1],TORSION_GROUP_REPS[t2],
                        TORSION_GROUP_REPS[t3],TORSION_GROUP_REPS[t4],
                        #ACLASS_ID[t1],ACLASS_ID[t2],
                        #ACLASS_ID[t3],ACLASS_ID[t4],
                        self.torsion_counts[i]))
    
mol2files = [l[:-1] for l in file(sys.argv[1])]

stat = StatClass()
classifier = AtomTypeClassifier()
torsassigner = TorsionAssigner()

stat.init_param(torsassigner.torsdb.db)

for mol2 in mol2files:
    molecule = MoleculeClass(mol2)
    classifier.apply_to_molecule(molecule)
    molecule.assign_rotable_torsions()
    torsassigner.assign(molecule,debug=False)
    stat.append_stat(molecule)

#stat.report_atypes()
#stat.report_torsions()
#stat.report_bond_params()
#stat.report_angle_params()
stat.report_torsion_params(torsassigner.torsdb.db)
