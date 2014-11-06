#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/Structure.py
## @brief  Class for structural representations of specific protein types.  Antibody and CDR work, feel free to add.  CDR Stuff Will be in C++ Rosetta soon.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
import os
import sys

class Antibody_Structure:
    """
    Simple class for accessing Modified_AHO antibody numbering information outside of Rosetta.  Use protocols/antibody/AntibodyInfo if importing Rosetta.
    """
    def __init__(self):
        #3 Different ways to access the CDR Loops.  Should be one.  Not convoluted.
        self.L1 = CDR("L1"); self.L2 = CDR("L2"); self.L3=CDR("L3");
        self.H1 = CDR("H1"); self.H2 = CDR("H2"); self.H3=CDR("H3");
        self.CDRS = [self.L1, self.L2, self.L3, self.H1, self.H2, self.H3]; #Switch to a dictionary with CDR name and cdr as item????
        

        self.CDR = {"L1":self.L1, "L2":self.L2, "L3":self.L3, "H1":self.H1, "H2":self.H2, "H3":self.H3}

    def get_CDRs(self):
        return self.CDRS
    
    def get_CDR(self, cdr_name):
        return self.CDR[cdr_name]
        
    def get_PDB_name(self, Path):
        root, pdb = os.path.split(Path)
        pdb = pdb.split('.')[0]
        self.pdb = pdb
        return self.pdb
    

class FRAMEWORK:
    def __init__(self):
        pass
    
    
class CDR:
    def __init__(self, name):
        self.regions = {
            "L1":['L', 24, 42],
            "L2":['L', 57, 72],
            "L3":['L', 107, 138],
            "H1":['H', 24, 42],
            "H2":['H', 57, 69],
            "H3":['H', 107, 138]}
        
        self.name = name
        self.region = self.regions[self.name]
        self.chain = self.region[0]
        self.Nter = self.region[1]
        self.Cter = self.region[2]
        self.residues = dict()
    
    def __str__(self):
        return str(self.regions[self.name])
    
    def get_pdb_chain(self):
        return self.regions[self.name][0]
        
    def get_pdb_start(self):
        return self.regions[self.name][1]
    
    def get_pdb_end(self):
        return self.regions[self.name][2]
        
    def set_gene(self, gene):
        self.gene = gene
        
    def add_residue(self, name, num):
        self.residues[num]=Residue(name, num)
        
class Residue:
    def __init__(self, name, num):
        self.name = name
        self.num = num
        self.createAtoms()
        self.createBonds()
    def createAtoms(self):
        pass
    def createBonds(self):
        pass
    
        #Should have Atoms, weight, it's own energy, the rotamer prob, etc.
    def caclulateDihedral(self, atom1, atom2, atom3):
        pass
############################################################FUTURE1#############################################################
class PDB:
    """
    The PDB should have seperate molecules for each chain.
    """
    def __init__(self, filename):
        self.loadPDB()
    def loadPDB(self):
        #Sets a few different properties.
        #Can look up residue data, atom data, etc.  =
        pass
        




#############################################################FUTURE2#############################################################    
class molecule:
    def __init__(self, name, atomsormolecule):
        pass
        #Were going to take all atoms one by one and create a molecule.
        #We can take two or more molecules and put them together to create a larger molecule.
        #Examples of a molecule:
            #H and L chain are each molecules.
            #Any ligands are molecules.
            #Molecules can covelently bond to other molecules.
            #Bonding can be a rule based on chemistry in the far flung future.
            

class protein:
    def __init__(self):
        pass
class protein_info:
    """
    The protein has protein molecules.  The PDB has a protein.  Need to write this carefully.
    """
    
    def __init__(self, sequence):
        #Were going to take the sequence.  Find the sequence's family.
        pass
        
    def giveStructure(self, pose):
        #Here, we give the sequence a structure.  We compare it to known structure and find it's family.  We parse it's components into domains using Pfam.
        #We get even more information on those domains. 
        pass
    
    def getInfo(self):
        #Here, we parse Uniprot.  We get all the information we possibly can.
        #We want to KNOW the species.  It's job.  It's fold.  It's required ligand.  Everything we can find on it.
        #Parsable.  Knowable.  
        pass
    
    def getPartners(self):
        #Here, we parse ProtCid. We return all known partners.  We try to determine whether the protein is a monomer or dimer, or something else.
        pass
    
    def breakIntoDomains(self):
        #Here, we break the structure into all domains.  This is useful if we want to do things to the domains.  OR start engineering them.
        pass
    def attachDomains(self, domain1, domain2):
        #Here we attach the two domains.  We will require more information about the design, but this will be added on later.
        pass
    
    
    
#If I ever need it:

class Atom:
    def __init__(self, name):
        pass
    def __add__(self, atom):
        """
        If we add two atoms together.  We create - A molecule, and a bond.
        This would be cool if we could get it to work.
        """
        return Bond(self, atom)
        #self.vdw
        #self.electrons
        #self.energy
        #self.
        pass
    def __iadd__(self):
        return Bond(self, Atom(self.name))
        
        
class Bond:
    def __init__(self, atom1, atom2):
        pass
    
    
    

#Just for Fucks sake:
class nucleus:
    pass
class electron:
    pass

    
