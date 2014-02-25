#!/usr/bin/env python

"""Make non-backbone heavy atom coordinate constraints from an
input PDB file. """

 
import os, sys

class Residue:
    def __init__(self,chain,resi,insert,resn):
        self.atoms_ = []
        self.chain_=chain
        self.resi_=resi
        self.insert_=insert
        self.resn_=resn
        self.pose_=None
    def compare(self,chain,resi,insert,resn):
        return (self.chain_==chain and self.resi_==resi and
                self.insert_==insert and self.resn_==resn)

    def add_atom(self, atom):
        self.atoms_.append(atom)
    def pose_num(self, num):
        self.pose_=num

class Atom:
    def __init__(self, name, x, y, z):
        self.name_ =name
        self.x_=x
        self.y_=y
        self.z_=z

def find_CA(residues):
    for res in residues:
        for atom in res.atoms_:
            if atom.name_.strip() == 'CA':
                return ' '.join((atom.name_, str(res.pose_)))
    return None
    
def main(filename, width, stdev):
    if filename.endswith('.pdb'):
        tag = filename[:-4]
    else:
        tag = filename

    chains = {}
    residues = []
    curres = None
    f = open(filename)
    try:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            atom = line[12:16]
            if atom[0] == 'H' or atom[1]=='H':
                continue # ignore hydrogens completely
            chain = line[21]
            resi = line[22:26]
            resn = line[17:20].upper()
            insert = line[26]
            x = line[30:38]
            y = line[38:46]
            z = line[46:54]
            if curres is None or not curres.compare(chain,resi,insert,resn):
                curres = Residue(chain,resi,insert,resn)
                chains.setdefault(chain,[]).append(curres)
                residues.append(curres)
            curres.add_atom(Atom(atom,x,y,z))
    finally:
        f.close()

    #Stupid hack to convert from PDB numbering to pose numbering
    chns = chains.keys()
    chns.sort()
    num = 1
    for c in chns:
        for i, r in sampling(chains[c]):
            r.pose_num(num)
            num += 1
            
    ref = None
    f = open(tag+'_sc.cst','w')
    try:
        for res in residues:
            resn = res.resn_
            pose = str(res.pose_)
            for atom in res.atoms_:
                sname = atom.name_.strip()
                name = atom.name_
                x = atom.x_
                y = atom.y_
                z = atom.z_
                
                if sname == 'CA':
                    ref = ' '.join((name, str(res.pose_)))
                if sname in ['C','CA','N','O','OXT']:
                    continue # ignore backbone atoms for constraints
                    # Note that 'O' also excludes water.
                if ref is None:
                    ref = find_CA(residues)
                    if ref is None:
                        raise ValueError("Pose must have CA *somwehere*!")
                
                if resn == 'ASN' and sname in ['OD1','ND2']: # Accomodate -flipHNQ
                    f.write("AmbiguousConstraint\n")
                    f.write(' '.join(["CoordinateConstraint",'OD1',pose,ref,x,y,z,
                                  "BOUNDED","0",width,stdev,"0.5","tag"])+'\n')
                    f.write(' '.join(["CoordinateConstraint",'ND2',pose,ref,x,y,z,
                                  "BOUNDED","0",width,stdev,"0.5","tag"])+'\n')
                    f.write("END_AMBIGUOUS\n")
                elif resn == 'GLN' and sname in ['OE1','NE2']: # Accomodate -flipHNQ
                    f.write("AmbiguousConstraint\n")
                    f.write(' '.join(["CoordinateConstraint",'OE1',pose,ref,x,y,z,
                                  "BOUNDED","0",width,stdev,"0.5","tag"])+'\n')
                    f.write(' '.join(["CoordinateConstraint",'NE2',pose,ref,x,y,z,
                                  "BOUNDED","0",width,stdev,"0.5","tag"])+'\n')
                    f.write("END_AMBIGUOUS\n")
                elif resn == 'HIS' and sname in ['ND1','CD2']: # Accomodate -flipHNQ
                    f.write("AmbiguousConstraint\n")
                    f.write(' '.join(["CoordinateConstraint",'ND1',pose,ref,x,y,z,
                                  "BOUNDED","0",width,stdev,"0.5","tag"])+'\n')
                    f.write(' '.join(["CoordinateConstraint",'CD2',pose,ref,x,y,z,
                                  "BOUNDED","0",width,stdev,"0.5","tag"])+'\n')
                    f.write("END_AMBIGUOUS\n")
                elif resn == 'HIS' and sname in ['CE1','NE2']: # Accomodate -flipHNQ
                    f.write("AmbiguousConstraint\n")
                    f.write(' '.join(["CoordinateConstraint",'CE1',pose,ref,x,y,z,
                                  "BOUNDED","0",width,stdev,"0.5","tag"])+'\n')
                    f.write(' '.join(["CoordinateConstraint",'NE2',pose,ref,x,y,z,
                                  "BOUNDED","0",width,stdev,"0.5","tag"])+'\n')
                    f.write("END_AMBIGUOUS\n")
                else:
                    f.write(' '.join(["CoordinateConstraint",name,pose,ref,x,y,z,
                                  "BOUNDED","0",width,stdev,"0.5","tag"])+'\n')
    finally:
        f.close()
        
            

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Usage: sidechain_cst.py <input.pdb> width stdev"
    else:
        main(sys.argv[1],sys.argv[2],sys.argv[3])
