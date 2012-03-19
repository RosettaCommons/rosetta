from Bio.PDB import *
from rosettautil.util import fileutil
class Atom:
    def __init__(self,name,type,chain,charge):
        self.name = name
        self.type = type
        self.chain = chain
        self.charge = charge
        self.bonds = []
    def add_bond(self,bond):
        self.bonds.append(bond)

class Bond:
    def __init__(self,lower,upper,type=0):
        self.lower = lower
        self.upper = upper
        self.type = type
    
class icoor:
    def __init__(self,name,phi,theta,d, stub1, stub2, stub3):
        self.name = name
        self.phi = phi
        self.theta = theta
        self.distance = d
        self.stub1 = stub1
        self.stub2 = stub2
        self.stub3 = stub3

class params:
    def __init__(self,filename):
        """load a param file into a data structure"""
        self.atoms = {}
        self.bonds = []
        self.icoors = {}
        infile = fileutil.universal_open(filename,'r')
        for line in infile:
            if len(line)==0:
                continue
            line = line.split()
            tag = line[0]
            if tag =="NAME":
                self.name = line[1]
            elif tag =="IO_STRING":
                self.io_string = line[1]
            elif tag == "TYPE":
                self.type = line[1]
            elif tag == "AA":
                self.aa = line[1]
            elif tag == "ATOM":
                current_atom = Atom(line[1],line[2],line[3],line[4])
                self.atoms[line[1]]= current_atom
            elif tag == "BOND":
                current_bond = Bond(line[1],line[2])
                self.bonds.append(current_bond)
                self.atoms[line[1]].add_bond(current_bond)
                self.atoms[line[2]].add_bond(current_bond)
            elif tag == "NBR_ATOM":
                self.nbr_atom = line[1]
            elif tag == "NBR_RADIUS":
                self.nbr_radius = line[1]
            elif tag == "ICOOR_INTERNAL":
                current_icoor = icoor(line[1],line[2],line[3],line[4],line[5],line[6],line[7])
                self.icoors[line[1]]=current_icoor
        infile.close()

    def get_atom(self,name):
        return self.atoms[name]

    def get_bond(self,lower,upper):
        for bond in self.bonds:
            if bond.lower == lower and bond.upper == upper:
                return bond

    def get_bonds_for_atom(self,name):
        atom = self.atoms[name]
        return atom.bonds

    def get_icoor(self,name):
        return self.icoors[name]
