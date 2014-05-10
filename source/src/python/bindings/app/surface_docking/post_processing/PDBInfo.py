#!/usr/bin/python
# Edited by Emily Koo to retain protein info and speed up processing
# Make sure to reset lists for each new PDB

from IO import load

class PDB:

    """
    Reduce a pdb atom or hetatm line to a set of objects
    describing that atom
    """

    def __init__(self, PDB):

        """
        Going column by column through pdb
        """
        
        #Atomic group i.e., atom or hetatm
        protein_or_surface = PDB[:6]
        if protein_or_surface == 'ATOM  ':
            self.atomic_group = 'protein'
        elif protein_or_surface == 'HETATM':
            self.atomic_group = 'surface'
          
        self.AtmNum = int(PDB[6:11]) 
        self.AtmTyp = str(PDB[11:16]).strip()
        self.Elem   = str(PDB[13:14])
        self.ResTyp = str(PDB[16:20]).strip()
        self.ResNum = int(PDB[22:26])
        self.X      = float(PDB[30:38])
        self.Y      = float(PDB[38:46])
        self.Z      = float(PDB[46:54])
        self.XYZ    = [float(PDB[30:38]),float(PDB[38:46]),float(PDB[46:54])]
#=================================================================

"""
Loads desired pdb info
"""
ATOMS = []
protein = []
surface = []
prot_coor = []
surf_coor = []

def resetlists():
 
    ATOMS[:] = []
    protein[:] = []
    surface[:] = []
    prot_coor[:] = []
    surf_coor[:] = []

def GetAtoms(pdb):

    if not ATOMS:
        for atom in pdb:
            if atom[:6] in ['ATOM  ', 'HETATM']:
                ATOM = PDB(atom)
                ATOMS.append(ATOM)
    return ATOMS

def partner(pdb):
    if not protein and not surface:
        ATOMS = GetAtoms(pdb)
        for ATOM in ATOMS:
            if ATOM.atomic_group == 'protein':
                protein.append(ATOM)
            if ATOM.atomic_group == 'surface':
                surface.append(ATOM)
                    
    return [protein, surface]

def coordinates(pdb):
    if not prot_coor or not surf_coor:
        ATOMS = GetAtoms(pdb)
        for ATOM in ATOMS:
            if ATOM.atomic_group == 'protein':
                coordinates = ATOM.X, ATOM.Y, ATOM.Z
                prot_coor.append(coordinates)
            if ATOM.atomic_group == 'surface':
                coordinates = ATOM.X, ATOM.Y, ATOM.Z
                surf_coor.append(coordinates)
                
    return [prot_coor, surf_coor]

def backbone(pdb):
    bbAtoms = ['N ', 'CA', 'C ', 'O ']
    backbone = []
    for atom in pdb:
        if atom[13:15] in bbAtoms:
            backbone.append(atom)

    return backbone

def GetResidueAtoms(Decoy, ResidueNumber):
    Residue = []

    protein = partner(load(Decoy))[0]
    for ATOM in protein:
        if ATOM.ResNum is ResidueNumber:
            Residue.append(ATOM)

    return Residue

def ResidueList(Decoy):
    ResidueList = []
    
    protein = partner(load(Decoy))[0]
    for ATOM in protein:
        if ATOM.ResNum not in ResidueList:
            ResidueList.append(ATOM.ResNum)

    return ResidueList

def GetAtomTypeInfo(Atoms,AtomType):
    MyAtoms = []

    for Atom in Atoms:
        if PDB(Atom).AtmTyp == AtomType:
            MyAtoms.append(Atom)

    return MyAtoms

def AminoAcidTable(ThreeLetterCode):
    OneLetterCode=''
    
    if ThreeLetterCode == 'ALA': OneLetterCode = 'A'
    if ThreeLetterCode == 'CYS': OneLetterCode = 'C'
    if ThreeLetterCode == 'ASP': OneLetterCode = 'D'
    if ThreeLetterCode == 'GLU': OneLetterCode = 'E'
    if ThreeLetterCode == 'PHE': OneLetterCode = 'F'
    if ThreeLetterCode == 'GLY': OneLetterCode = 'G'
    if ThreeLetterCode == 'HIS': OneLetterCode = 'H'
    if ThreeLetterCode == 'ILE': OneLetterCode = 'I'
    if ThreeLetterCode == 'LYS': OneLetterCode = 'K'
    if ThreeLetterCode == 'LEU': OneLetterCode = 'L'
    if ThreeLetterCode == 'MET': OneLetterCode = 'M'
    if ThreeLetterCode == 'ASN': OneLetterCode = 'N'
    if ThreeLetterCode == 'SEP': OneLetterCode = 'O'
    if ThreeLetterCode == 'PRO': OneLetterCode = 'P'
    if ThreeLetterCode == 'GLN': OneLetterCode = 'Q'
    if ThreeLetterCode == 'ARG': OneLetterCode = 'R'
    if ThreeLetterCode == 'SER': OneLetterCode = 'S'
    if ThreeLetterCode == 'THR': OneLetterCode = 'T'
    if ThreeLetterCode == 'VAL': OneLetterCode = 'V'
    if ThreeLetterCode == 'TRP': OneLetterCode = 'W'
    if ThreeLetterCode == 'TYR': OneLetterCode = 'Y'
    if ThreeLetterCode == 'CGU': OneLetterCode = 'Z'

    return OneLetterCode


