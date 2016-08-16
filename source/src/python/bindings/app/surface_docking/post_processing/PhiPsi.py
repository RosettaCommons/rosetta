#!/usr/bin/env python

from string import strip, split
from sys import argv, exit
from math import sqrt, pi, acos, cos, sin
#from numpy import zeros, arange

class ATOMinfo:

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


        #Atomic coordinates
        self.X = float(PDB[30:38])
        self.Y = float(PDB[38:46])
        self.Z = float(PDB[46:54])

#=================================================================

"""
Loads desired pdb info
"""
def load(file):
    file = open(file, 'r').readlines()
    loaded_file = []
    for line in file:
        loaded_file.append(strip(line))

    return loaded_file

def partner(pdb):
    protein = []
    surface = []
    for atom in pdb:
        if atom[:6] in ['ATOM  ', 'HETATM']:
            ATOM = ATOMinfo(atom)
            if ATOM.atomic_group == 'protein':
                protein.append(atom)
            if ATOM.atomic_group == 'surface':
                surface.append(atom)

    return [protein, surface]

def coordinates(pdb):
    protein = []
    surface = []
    for atom in pdb:
        if atom[:6] in ['ATOM  ', 'HETATM']:
            ATOM = ATOMinfo(atom)
            if ATOM.atomic_group == 'protein':
                coordinates = ATOM.X, ATOM.Y, ATOM.Z
                protein.append(coordinates)
            if ATOM.atomic_group == 'surface':
                coordinates = ATOM.X, ATOM.Y, ATOM.Z
                surface.append(coordinates)

    return [protein, surface]

def backbone(pdb):
    bbAtoms = ['N ', 'CA', 'C ']
    backbone = []
    for atom in pdb:
        if atom[13:15] in bbAtoms:
            backbone.append(atom)

    return backbone

def CA(pdb):
    CAAtoms = ['CA']
    CA = []
    for atom in pdb:
        if atom[13:15] in CAAtoms:
            CA.append(atom)

    return CA

def hbond_record(pdb):
    donor = []
    acceptor = []
    Pair = []
    Type = []
    interface_protein_atom = []
    end = 1
    appending_protein_record = 0
    for line in pdb:
	if len(split(line)) == 14 and split(line)[0] != '#Dch':

             donor.append(split(line)[1:4]) # Get columns of index start:end+1
             #print "donor ", split(line)[1:4]
             acceptor.append(split(line)[5:8])
             #print "acceptor ",split(line)[5:8]

	"""
        if not end or line[:5] == 'begin':
	   #print split(line)[0]
            if line[:4] == 'end ':
                end = 1
                appending_protein_record = 0
                appending_surface_record = 0
            elif appending_protein_record or split(line)[1] == 'protein':
                end = 0
                appending_protein_record = 1
                if split(line)[0] != 'begin' and split(line)[0] != 'donor' and split(line)[0] != '#Dch' and split(line)[0] != ' ':
		    if len(split(line)) == 14:
	                    donor.append(split(line)[1:4])
			   #print "donor ", split(line)[1:4]
	                    acceptor.append(split(line)[5:8])
			   #print "acceptor ",split(line)[5:8]
          """

    for residue in range(len(donor)):
        type = acceptor[residue][2], donor[residue][2]
        Type.append(type)
        pair = int(acceptor[residue][0]), int(donor[residue][0])
        Pair.append(pair)

    return Pair, Type

#======================================================================

"""
A bunch of vector equations
"""
def vector(pointA, pointB):
    vector = [0,0,0]
    vector[0] = pointB[0] - pointA[0]
    vector[1] = pointB[1] - pointA[1]
    vector[2] = pointB[2] - pointA[2]

    return vector

def length(vector):
    length = sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)

    return length

def unit_vector(vector, length):
    unit_vector = [vector[0]/length, vector[1]/length, vector[2]/length]

    return unit_vector

def dot_product(vector1, vector2):
    dot_product = vector1[0]*vector2[0] + vector1[1]*vector2[1] + vector1[2]*vector2[2]

    return dot_product

def cross_product(vector1, vector2):
    cross_product = [0,0,0]
    cross_product[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1]
    cross_product[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2]
    cross_product[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0]

    return cross_product

def angle(vector1, vector2):
    unit_vector1 = unit_vector(vector1, length(vector1))
    unit_vector2 = unit_vector(vector2, length(vector2))
    angle = acos(dot_product(unit_vector1, unit_vector2))*180/pi

    return angle

def scalar_triple_product(vector1, vector2, vector3):
    scalar_triple_product = dot_product( vector1, cross_product(vector2, vector3) )

    return scalar_triple_product

#============ Get The Info =============================

def PhiPsi(pdb):
    Phi = []
    Psi = []
    protein_backbone = backbone(partner(load(pdb))[0])
    protein_backbone_coordinates = coordinates(backbone(partner(load(pdb))[0]))[0]
    Phi.append(360.0)
    for atom in range(len(protein_backbone) - 1):
        if protein_backbone[atom][13:15] in ['C ']:
            phi = angle(
            cross_product(
                vector(protein_backbone_coordinates[atom], protein_backbone_coordinates[atom + 1]), vector(protein_backbone_coordinates[atom + 1], protein_backbone_coordinates[atom + 2])),
            cross_product(
                vector(protein_backbone_coordinates[atom + 1], protein_backbone_coordinates[atom + 2]), vector(protein_backbone_coordinates[atom + 2], protein_backbone_coordinates[atom + 3])))
            if ( scalar_triple_product( vector(protein_backbone_coordinates[atom + 2], protein_backbone_coordinates[atom + 3]),
                                        vector(protein_backbone_coordinates[atom], protein_backbone_coordinates[atom + 1]),
                                        vector(protein_backbone_coordinates[atom + 1], protein_backbone_coordinates[atom + 2]))
                < 0.0 ):
                phi = -phi
            else:
                phi = phi
	   #print "phi ", phi
            Phi.append(phi)

    for atom in range(len(protein_backbone) - 3):
        if protein_backbone[atom][13:15] in ['N ']:
            psi = angle(
            cross_product(
                vector(protein_backbone_coordinates[atom], protein_backbone_coordinates[atom + 1]), vector(protein_backbone_coordinates[atom + 1], protein_backbone_coordinates[atom + 2])),
            cross_product(
                vector(protein_backbone_coordinates[atom + 1], protein_backbone_coordinates[atom + 2]), vector(protein_backbone_coordinates[atom + 2], protein_backbone_coordinates[atom + 3])))
            if ( scalar_triple_product( vector(protein_backbone_coordinates[atom + 2], protein_backbone_coordinates[atom + 3]),
                                        vector(protein_backbone_coordinates[atom], protein_backbone_coordinates[atom + 1]),
                                        vector(protein_backbone_coordinates[atom + 1], protein_backbone_coordinates[atom + 2]))
                < 0.0 ):
                psi = -psi
            else:
                psi = psi
	   #print "psi ", psi
            Psi.append(psi)

    Psi.append(360.0)
    return [Phi, Psi]

def SecondaryStructure(pdb):
    Phi = PhiPsi(pdb)[0]
    Psi = PhiPsi(pdb)[1]
    Helix = []
    Loop = []
    SecStruct = []
    HBond = hbond_record(load(pdb))[0]
    HBondType = hbond_record(load(pdb))[1]
    Protein = CA(partner(load(pdb))[0])
    AllowedHBondDist_Helix = 4
    AllowedHBondDist_Turns = [3,4,5]

    '''============= Detect Helical Residues ==============='''

    HelixBin = []
    for hbond in range(len(HBond)):
        acceptorRes = HBond[hbond][0]
        donorRes = HBond[hbond][1]
        acceptorType = HBondType[hbond][0]
        donorType = HBondType[hbond][1]
        HBondDist = abs(donorRes - acceptorRes)
        if acceptorType == 'O' and donorType == 'N' and HBondDist == AllowedHBondDist_Helix:
            HelixBin.append(HBond[hbond][0])
    HelicalSegmentBin = []
    for Residue in range(len(Protein)):
        if int(Protein[Residue][24:26]) in HelixBin and int(Protein[Residue][24:26]) - 1 in HelixBin:
            for HelixResidue in range(int(Protein[Residue][24:26]), int(Protein[Residue][24:26]) + 4):
                HelicalSegmentBin.append(HelixResidue)
    for Residue in Protein:
        if int(Residue[24:26]) in HelicalSegmentBin:
            Helix.append(1)
        else:
            Helix.append(0)



      #HelixVar=25.0
    #for residue in range(len(Protein)):
    #    if residue != len(Protein) - 1:
    #        if int(Psi[residue] + Phi[residue + 1]) in range(-130, -80):
    #            Helix.append(1)
    #        else:
    #            Helix.append(0)
    #    else:
    #        Helix.append(0)
       #if Phi[residue] > -64.0 - HelixVar and Phi[residue] < -64.0 + HelixVar and Psi[residue] > -41.0 - HelixVar and Psi[residue] < -41.0 + HelixVar:
    #    Helix.append(1)
        #else:
        #    Helix.append(0)

    '''============= Detcect Loop Residues ==============='''
    LoopBin = []
    for hbond in range(len(HBond)):
        acceptorRes = HBond[hbond][0]
        donorRes = HBond[hbond][1]
        HBondDist = abs(donorRes - acceptorRes)
        startRes = min(donorRes, acceptorRes)
        if HBondDist in AllowedHBondDist_Turns:
            for Residue in range(startRes, startRes + HBondDist + 1):
                LoopBin.append(Residue)
                #HBondMatrix[donor, Residue] += 1

    #ndex = [arange(0,len(Protein),1),arange(0,len(Protein),1)]

    #for x in enumerate(Index[0]):
    #    for y in enumerate(Index[1]):
    #        if HBondMatrix[x[1], y[1]] != 0:
    #           #print x[1],' ', y[1], ' ', HBondMatrix[x[1],y[1]]

    for Residue in Protein:
        if int(Residue[24:26]) in LoopBin:
            Loop.append(1)
        else:
            Loop.append(0)

    #while residue in range(len(Protein)):
    #    ResIsHBonded=0
    #    TmpBin = []
    #    ResNum = residue + 1
    #    for hbond in range(len(HBond)):
    #        if ResNum == HBond[hbond][0] or ResNum ==  HBond[hbond][1]:
    #            HBondSpace = abs(HBond[hbond][1] - HBond[hbond][0])
    #            if HBondSpace in AllowedPrimaryDistanceBetweenHBonds:
    #                TmpBin.append(HBondSpace)
    #                ResIsHBonded=1


            #if hbond == len(HBond) -1:
     #   if ResIsHBonded:
     #       HBondSpace = max(TmpBin)
     #       for loop_residue in range(HBondSpace + 1):
     #           residue += 1
     #           Loop.append(1)

     #   else:
     #       residue += 1
     #       Loop.append(0)

    #print len(Loop)
    '''======= Determine Loop, Helical, or Extended segments'''

    residue=0
    while residue in range(len(Protein)):

        '''========== First Helical =============='''

        if Helix[residue]:
            SecStruct.append('Helix')
            #residue += 1
      #  helical_segment=0
      #  if Helix[residue]:
      #      if residue!= 0 and SecStruct[residue - 1] == 'Helix':
      #          SecStruct.append('Helix')
      #          helical_segment=1
      #          residue += 1
      #      else:
      #          for helical_residue in range(residue, residue + 5):
      #              if not Helix[helical_residue]:
      #                  helical_segment=0
      #                  break
      #              else:
      #                  helical_segment=1
      #                  continue
      #          if helical_segment:
      #              for helical_residue in range(residue, residue + 5):
      #                  SecStruct.append('Helix')
      #                  residue += 1

        '''========== Then Loop ==================='''

        if residue < len(Protein):
            if not Helix[residue] and Loop[residue]:
            #if not helical_segment and Loop[residue]:
                SecStruct.append('Loop')
                #residue += 1

        '''======= Then Extended  ================='''

        if residue < len(Protein):
            if not Helix[residue] and not Loop[residue]:
            #if not helical_segment and not Loop[residue]:
                SecStruct.append('Extended')
                #residue += 1
        residue += 1
    return SecStruct

#=========== Make The Output =====================

def MAIN(pdb):
    Phi = PhiPsi(pdb)[0]
    Psi = PhiPsi(pdb)[1]
    SecStruct = SecondaryStructure(pdb)
    FileName = argv[1][:-4] + '.dss'
    OpenFile = open(FileName, 'w')
    Protein = CA(partner(load(pdb))[0])
    OpenFile.write("Resn Resi     Phi     Psi    SS Type" + '\n')
    for residue in range(len(Protein)):
        Line = "%4s %4s %7.6s %7.6s %1s %1s"%(Protein[residue][22:26], Protein[residue][17:20], Phi[residue], Psi[residue], '', str.ljust(SecStruct[residue], 1))
        OpenFile.write(Line + '\n')


pdb = argv[1]
MAIN(pdb)

