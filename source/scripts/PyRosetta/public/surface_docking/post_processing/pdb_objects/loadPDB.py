#!/usr/bin/env python

##############################################
#Module loadPDB.py
#created 10/13/03, Mike Daily
#last modified 3/22/04, Mike Daily
##############################################
import string
import coordlib
import sys
import copy
import os
#from numarray import array, matrixmultiply
from numpy import numarray
from math import sin, cos, sqrt

#*****************************************************************************
#FUNCTIONS FOR LOADING PDB INTO PDB CLASS

output_clean_pdb = 0

def corr_int(resid):
    resid_num = ''
    for char in resid:
        if char in string.digits:
            resid_num = resid_num + char
    return int(resid_num)

class PDBatom:
    """
    reduce a pdb atom line (columns of text) to a set of objects
    describing that atom
    """
    def __init__(self, PDBline):
        recordType=PDBline[0:4]
        self.atomtype = string.strip(PDBline[13:16])
        self.restype = string.strip(PDBline[17:20])
        self.chain = PDBline[21]
        #Add in a dummy chain if there is none from the pdb file
        #(Makes later data processing easier)
        if self.chain == ' ':
            self.chain = '0'
        self.resid = string.strip(PDBline[22:27])
        #Please note resid is a STRING not an int (works better
        #if resid is something like '14A')
        self.coords = [float(PDBline[30:38]), float(PDBline[38:46]),
                       float(PDBline[46:54])]
        if recordType == 'ATOM':
            #recognize DNA residues
            if self.restype in 'AGCT':
                self.groupType = 'dna'
            else:
                self.groupType='std'
        elif recordType == 'HETA':
            self.groupType='het'


        #temperature factor
        try:
            self.Bfac = float(string.strip(PDBline[60:65]))
        except ValueError:
            self.Bfac = 0.0

class new_PDBatom:
    #For creation of PDBatoms from coordinates
    def __init__(self, groupType, atomtype, restype, chain, resid, coords):
        self.groupType = groupType
        self.atomtype = atomtype
        self.restype = restype
        self.chain = chain
        self.resid = resid
        self.coords = coords

#============================================================================

class PDBres:
    """
    take a resatomlist (output of getResidues) and reduces it to a
    residue class, consisting of a type, chain, resid, and a list of
    atoms.  Each atom has been reduced from a pdbatom (output of class
    PDBatom) to a reduced atom (output of reducedAtom).  The pdbatom
    information is no longer needed, since it is the same for every residue.

    All the methods below are part of this class, down to the next line of
    '=' characters.
    """
    def __init__(self, resatomlist):
        self.groupType=resatomlist[0].groupType
        self.restype = resatomlist[0].restype
        try:
            self.resaa = restypedict[self.restype]
        except KeyError:
            self.resaa = 'NA'
        self.chain = resatomlist[0].chain
        self.resid = resatomlist[0].resid
        self.resnum = corr_int(self.resid)
        self.setPrintId()

        #name can be reset if I want to change the name of a residue for
        #computing purposes but not its output label.
        self.name = self.printid
        self.Bfac = 0.0
        self.initializeAtoms(resatomlist)
        self.max_radius = None
        self.med_radius = None

    def setPrintId(self):
        self.printid = self.restype + ' ' + self.chain + ' ' + self.resid

    def setChain(self, chain):
        self.chain = chain
        self.setPrintId()

    def initializeAtoms(self, resatomList):
        """
        Initializes several definitions of the atoms in resatomlist:

        1) direct-callable:  calling on the name of the atom records
        coordinates of that atom
        2) list of atoms:  list of coordinates; for blind iteration over
        atom coordinates (e.g. calculating centroid)
        3) list of types:  which atom types does the residue contain
        """
        self.coordlist = []
        self.atomlist = []
        for atom in resatomList:
            #don't load in multiple copies of the same atom (e.g. CB A, CB A)
            if not(atom.atomtype in self.atomlist):
                #initializes objects for direct calling of atoms (e.g. res.N = ?)
                coords = atom.coords
                setattr(self, atom.atomtype, coords)
                #makes coordlist, atomtypelist, namedAtomlist
                self.coordlist.append(coords)
                self.atomlist.append(atom.atomtype)

                #Get B-factor from CA atom
                if atom.atomtype == 'CA':
                    self.Bfac = atom.Bfac
        self.Natoms = len(self.atomlist)

    #I have made getitem equivalent to getattr so that I can easily reference
    #atoms using the [] overloading convention, e.g. res['N'] if I want to
    #loop over a list of atomtypes without using getattr(res, atomtype)

    def hasatom(self, atomtype):
        if atomtype in self.atomlist:
            return 1
        return 0

    def copy(self):
        tmp_copy = copy.copy(self)
        tmp_copy.coordlist = self.coordlist[:]
        tmp_copy.atomlist = self.atomlist[:]
        return tmp_copy

    def __getitem__(self, name):
        return getattr(self, name)

    def __setitem__(self, name, value):
        return setattr(self, name, value)

    global bbAtomList
    bbAtomList = ['N', 'CA', 'C', 'O']

    def getScAtoms(self):
        scAtomList = []
        for atomtype in self.atomlist:
            if not (atomtype in bbAtomList):
                scAtomList.append(atomtype)
        self.scatomlist = scAtomList

    def addCentroid(self):
        #adds centroid and max radius (often needed for calc. with centroid)
        self.centroid = coordlib.vavg(self.coordlist)
        cen_radii = []
        for atom in self.coordlist:
            cen_radii.append(coordlib.dist(atom, self.centroid))
        cen_radii.sort()
        self.max_radius = cen_radii[-1]

        #median radius:  a radius from the centroid enclosing half the atoms
        #of the residue
        middle = self.Natoms / 2
        if self.Natoms % 2:
            self.med_radius = cen_radii[middle]
        else:
            self.med_radius = (cen_radii[middle-1] + cen_radii[middle]) / 2.0

    def fetchAtoms(self, desiredAtomList):
        """
        input:  list of desired atom names from this residue
        output:  list of coordinates in the order of the list of
        desired atom names.
        If one of the atoms is not found, an error message is printed and
        an empty list is returned.
        The program which called this function can then interpret the empty
        list,which could mean a few things:
        1) coding error:  searched for atom types inappropriate for this
        residue type
        2) incomplete sidechain or backbone
        """
        coordlist = []
        for atomname in desiredAtomList:
            try:
                #get the coordinates
                coordlist.append(self[atomname])
            #atom not found in this residue!
            except AttributeError:
                sys.stderr.write('atom ' + atomname + ' not found in \
                residue ' + self.printid + '\n' )
                return []
        return coordlist

    def fetchBB(self):
        return self.fetchAtoms(bbAtomList)

    def fetchSC(self):
        if not(hasattr(self, 'scatomlist')):
            self.getScAtoms()
        return self.fetchAtoms(self.scatomlist)

    def addscCentroid(self):
        scatoms = self.fetchSC()
        #For glycines, set the sidechain centroid to the whole residue
        #centroid
        if scatoms == []:
            self.addCentroid()
            self.sccentroid = self.centroid
        else:
            self.sccentroid = coordlib.vavg(scatoms)
#==========================================================================
#Protein class for breaking a protein down into chains and residues

class ProteinChain:
    def __init__(self):
        self.pdbstart = None
        self.polypep = None
        self.ligands = None
        self.chain = None
        self.Nligands = None
        self.sequence = None
        self.unmapped_residues = []
        self.type = 'het'

    def parse(self, pdbstart):
        #The main goal here is to separate protein residues from ligand
        #residues.
        self.pdbstart = pdbstart
        self.polypep = []
        self.ligands = []
        self.chain = self.pdbstart[0].chain

        for res in self.pdbstart:
            if res.groupType in ['std', 'dna']:
                self.polypep.append(res)
            else:
                self.ligands.append(res)

        if self.polypep != None:
            self.type = self.polypep[0].groupType

        self.Nligands = len(self.ligands)
        self.sequence = getSequence(self.polypep)

    def __len__(self):
        return len(self.polypep)

    def __getitem__(self, i):
        #allows iteration over protein residues while ignoring ligands if
        #you want to
        return self.polypep[i]

    def apply_map(self, map):
        #map is a list of 1's and zero's corresponding to the residues in loadedPDB.
        #Keep residues corresponding to 1's and discard residues corresponding to
        #0's (the original pdb is lost)
        #print map
        newPDB = []
        Nres = len(self.polypep)
        #Make sure map is same length as pdb
        if Nres != len(map):
            print Nres, 'residues; map length=', len(map)
            print 'pdb length and map length do not match; exiting'
            sys.exit(1)

        lowest_true = 0
        highest_true = 0
        for i in range(Nres):
            if map[i] == 1:
                newPDB.append(self.polypep[i])
                if lowest_true == 0:  lowest_true = i
                highest_true = i

        #Iterate again to get unmapped residues
        #I only want internally non-matching residues, not non-matching residues
        #at the ends (may be an artifact of crystallography)
        #That is why I recorded the lowest mapping residue and highest mapping
        #residue in the previous loop
        self.unmapped_residues = []
        """
        for i in range(Nres):
            if map[i] == 0 and i > lowest_true and i < highest_true:
                self.unmapped_residues.append(self.polypep[i])
        """
        for i in range(Nres):
            if map[i] == 0:
                self.unmapped_residues.append(self.polypep[i])
        #replace self.polypep with the new, aligned PDB.
        self.polypep[:] = newPDB

    def transformall(self, tran, rot):
        self.polypep[:] = transformCoords(self.polypep, tran, rot)
        self.ligands[:] = transformCoords(self.ligands, tran, rot)

    def addCentroids(self):
        for res in self.polypep:
            res.addCentroid()
        for res in self.ligands:
            res.addCentroid()

    def copy(self):
        copy_chain = ProteinChain()
        copy_chain.polypep = self.polypep[:]
        copy_chain.ligands = self.ligands[:]
        copy_chain.chain = self.chain
        copy_chain.Nligands = self.Nligands
        copy_chain.sequence = self.sequence
        copy_chain.unmapped_residues = self.unmapped_residues
        return copy_chain

class chainType:
    def __init__(self, sequence, refChain, chainList):
        #sequence:  sequence of that chain type
        #refChain:  1st chainID in the protein having that sequence
        #chainList:  all chainIDs in the protein with that sequence
        self.sequence = sequence
        self.refChain = refChain
        self.chainList = chainList
        self.chainLabel = ''

    def setChainLabel(self):
        for chainID in self.chainList:
            self.chainLabel += chainID

class Protein:
    def __init__(self):
        self.pdbstart = None
        self.chainList = None
        self.chains = None
        self.unassigned_ligands = None
        self.Nchains = None
        self.typeList = None
        self.chainTypeList = None
        self.Ntypes = None
        self.unmapped_residues = []
        self.DNAchains = []

    def parse(self, pdbstart):
        DNAresidues = []
        self.pdbstart = []
        for res in pdbstart:
            if res.groupType == 'dna':
                DNAresidues.append(res)
            else:
                self.pdbstart.append(res)

        if len(DNAresidues) > 0:
            DNApdbs = parseByChain(DNAresidues)
            for pdb in DNApdbs:
                new_chain = ProteinChain()
                new_chain.parse(pdb)
                self.DNAchains.append(new_chain)

        #self.pdbstart = pdbstart
        self.chainList = [self.pdbstart[0].chain]
        self.chains = []
        self.unassigned_ligands = []

        Nres = len(self.pdbstart)
        last_chain = self.pdbstart[0].chain
        current_chain_resList = [self.pdbstart[0]]

        for i in range(1, Nres):
            current_res = self.pdbstart[i]
            current_chain = current_res.chain
            #Separate out unassigned ligands -
            #Otherwise they will mess up the assignment of chains
            if last_chain != '0' and current_chain == '0':
                self.unassigned_ligands.append(current_res)
                continue

            if current_chain == last_chain:
                current_chain_resList.append(current_res)
            else:
                new_chain = ProteinChain()
                new_chain.parse(current_chain_resList)
                self.chains.append(new_chain)
                self.chainList.append(current_chain)
                current_chain_resList = [current_res]
            last_chain = current_chain

        #append what's left over at the end
        new_chain = ProteinChain()
        new_chain.parse(current_chain_resList)
        self.chains.append(new_chain)
        self.Nchains = len(self.chainList)

        self.typeList = [chainType(self.chains[0].sequence, self.chainList[0],
                         [self.chainList[0]])]
        self.chainTypeList = [0] * self.Nchains

        for i in range(1, self.Nchains):
            Ntypes = len(self.typeList)
            matched = 0
            for j in range(Ntypes):
                #if sequence matches a known sequence
                if self.chains[i].sequence == self.typeList[j].sequence:
                    matched = 1
                    self.typeList[j].chainList.append(self.chainList[i])
                    self.chainTypeList[i] = j
                    break

            if matched:
                continue

            #create new type
            new_type = chainType(self.chains[i].sequence, self.chainList[i],
                                 [self.chainList[i]])
            self.typeList.append(new_type)
            #Make note of new type in chaintypelist at that position
            self.chainTypeList[i] = Ntypes

        for type in self.typeList:
            type.setChainLabel()

        self.Ntypes = len(self.typeList)
        #for type in self.typeList:
        #    print type.sequence, type.refChain, type.chainList, len(type.sequence)

        """
        #Condense types together with identical sequences except for gaps and
        #termini differences
        if self.Ntypes > 1:
            newChainTypeList = [self.typeList[0]]
            for i in range(1, self.Ntypes):
                aligned_sequences = alignPDB()
                aligned_sequences.
        """

    def getChain(self, chainID):
        #print self.chainList
        for i in range(self.Nchains):
            if self.chainList[i] == chainID:
                return self.chains[i]
        print 'chain', chainID, 'not found; exiting'
        sys.exit(1)

    def getAllProteinRes(self):
        """
        Returns a single list containing all protein (but not ligand) residues
        of the protein
        """
        wholeProtein = []
        for chain in self.chains:
            wholeProtein = wholeProtein + chain.polypep
        return wholeProtein

    def getAllRes(self):
        """
        Returns a single list containing all protein and ligand residues
        of the protein
        """
        wholeProtein = []
        for chain in self.chains:
            wholeProtein = wholeProtein + chain.polypep + chain.ligands
        wholeProtein = wholeProtein + self.unassigned_ligands
        return wholeProtein

    def getLigands(self):
        ligands = []
        for chain in self.chains:
            ligands = ligands + chain.ligands
        ligands = ligands + self.unassigned_ligands
        return ligands

    def removeLigands(self):
        for chain in self.chains:
            chain.ligands = []
        self.unassigned_ligands = []

    def apply_map(self, chain, map):
        #apply a map to a given chain of the protein
        #(keep residues with 1's and skip residues with zero's)
        #find the chain
        for i in range(self.Nchains):
            if self.chainList[i] == chain:
                #apply map to that chain
                self.chains[i].apply_map(map)
                self.unmapped_residues += self.chains[i].unmapped_residues
                return

    def transformall(self, tran, rot):
        newProtein = self.copy()
        for chain in newProtein.chains:
            chain.transformall(tran, rot)
        newProtein.unassigned_ligands[:] = transformCoords(
            newProtein.unassigned_ligands, tran, rot)
        return newProtein

    def addCentroids(self):
        for chain in self.chains:
            chain.addCentroids()
        for chain in self.DNAchains:
            chain.addCentroids()
        for res in self.unassigned_ligands:
            res.addCentroid()

    def copy(self):
        copy_protein = Protein()
        copy_protein.chainList = self.chainList[:]
        copy_protein.unassigned_ligands = self.unassigned_ligands[:]

        copy_protein.chains = []
        for chain in self.chains:
            copy_protein.chains.append(chain.copy())

        copy_protein.Nchains = self.Nchains
        copy_protein.typeList = self.typeList[:]
        copy_protein.chainTypeList = self.chainTypeList[:]
        copy_protein.Ntypes = self.Ntypes
        return copy_protein


#==========================================================================

restypedict = {
    'ALA':  'A',
    'CYS':  'C',
    'ASP':  'D',
    'GLU':  'E',
    'PHE':  'F',
    'GLY':  'G',
    'HIS':  'H',
    'ILE':  'I',
    'LYS':  'K',
    'LEU':  'L',
    'MET':  'M',
    'ASN':  'N',
    'PRO':  'P',
    'GLN':  'Q',
    'ARG':  'R',
    'SER':  'S',
    'THR':  'T',
    'VAL':  'V',
    'TRP':  'W',
    'TYR':  'Y',
    }

def is_bb_complete(PDBres):
    #Crashes the script if a residue in the PDB has a missing backbone
    atomlist=PDBres.atomlist
    bb_complete = (('N' in atomlist) * ('CA' in atomlist) *
                   ('C' in atomlist) * ('O' in atomlist))
    if not bb_complete:
        print 'loadPDB.py has found a missing backbone in PDB residue', \
              PDBres.printid
        print 'generating a new pdb without the missing backbone residue'
        print '(same name as original pdb)'
        print 'The old pdb has been moved to <pdb>.badbb.pdb'
    return bb_complete

def loadAtoms(pdb):
    """
    takes opened pdb file 'pdb' and reduces each ATOM entry to a
    PDBatom as defined
    in the PDBatom class above.
    the output is a list of PDBatoms
    At this step, waters and hydrogens are removed
    """
    atomlist = []
    for item in pdb:
        if item[0:4] in ['ATOM', 'HETA']:
            newatom = PDBatom(item)
            #remove waters and hydrogens!
            if newatom.restype != 'HOH' and newatom.atomtype[0] != 'H':
                atomlist.append(newatom)
    return atomlist

def parseAtomList(pdbatoms):
    """
    takes a list of pdbatoms (output of loadAtoms above) and splits
    it into sublists, one for each residue
    The output looks like the following:
    [[atom1, atom2, ... , atom-Natoms(res1)] , [atom1, atom2, ... ,
    atom-Natoms(res2)] , ... , Nres]
    The format of the pdbatom entries does not change, only the list
    structure.
    """
    if len(pdbatoms) == 0:
        return []
    outputlist = []
    resatomlist = [pdbatoms[0]]
    for i in range (1, len(pdbatoms)):
        current = pdbatoms[i]
        previous = pdbatoms[i-1]
        if [current.chain, current.resid] == [previous.chain,
                                              previous.resid]:
            resatomlist.append(current)
        else:
            outputlist.append(resatomlist)
            resatomlist = [current]
    outputlist.append(resatomlist)
    return outputlist

def convert_to_loadedPDB(resAtomList):
    """
    Input:  a list of PDBatoms parsed into residues
    Output: a loadedPDB (each residue is converted to a PDBres class instance)

    Residues with missing backbone atoms are skipped, and the global flag
    output_clean_pdb is turned on so that a new pdb will be created without the
    missing backbone
    """
    global output_clean_pdb #Set to 0 at beginning of file
    PDBresList = []
    for res in resAtomList:
        newRes = PDBres(res)
        if newRes.groupType == 'std' and not(is_bb_complete(newRes)):
            output_clean_pdb = 1
        else:
            PDBresList.append(newRes)
    return PDBresList

def loadProt(pdbfile):
    """
    combines above functions
    Input is a pdb file
    Output is a the output of function reduceReslist above.
    Function proceeds in several steps:
    1) Extract ATOM lines and reduce them to pdbatom classes (function
    loadAtoms)
    2) Parse pdbatom list (result of (1)) into residues (function
    getResidues) and convert to PDBres classes
    """
    global output_clean_pdb
    pdb=open(pdbfile, 'r').readlines()
    pdbid = string.split(pdbfile, '.')[0]

    pdbatomlist = loadAtoms(pdb)
    resAtomList = parseAtomList(pdbatomlist)
    loadedPDB = convert_to_loadedPDB(resAtomList)
    addIndices(loadedPDB)
    if output_clean_pdb:
        bad_pdb_filename = pdbid + '.badbb.pdb'
        os.system('mv ' + pdbfile + ' ' + bad_pdb_filename)
        writeSinglePDB(loadedPDB, pdbfile)
        output_clean_pdb=0
    return loadedPDB

def loadProtein(pdbfile):
    #loads a protein from a pdb into a loadProt object, then converts to
    #protein format (split into chains, heteroatoms, etc)
    loadedPDB = loadProt(pdbfile)
    newProt = Protein()
    newProt.parse(loadedPDB)
    return newProt

#*******************************************************************************
#FUNCTIONS FOR OUTPUTTING LOADEDPDB TO A FILE

def outputCoords(coords):
    #Converts coordinates from list [x,y,z] to string (PDB) format
    outstr = ''
    for item in coords:
        outstr = outstr + string.rjust(('%4.3f' % (item)),8)
    return outstr

def generatePDB(loadedPDB):
    """
    Output a loadedPDB in PDB format
    """
    outputList = []
    atomCounter = 1
    for res in loadedPDB:
        for atom in res.atomlist:
            recordType = 'ATOM'
            if res.groupType == 'het':
                recordType = 'HETATM'
            pdbline = string.ljust(recordType, 6) + \
                      string.rjust(str(atomCounter),5) + '  ' + \
                      string.ljust(atom, 4) + string.rjust(res.restype, 3) + ' ' + res.chain + \
                      string.rjust(str(res.resnum),4) + 4*' ' + \
                      outputCoords(res[atom]) + '  1.00\n'
            outputList.append(pdbline)
            atomCounter = atomCounter + 1
    return outputList

def generate_one_model(splitPDB):
    """
    Input:  a loadedPDB split by chains
    Output:  a list of PDB lines, with secondary structure at the beginning
    and TERs at the end of each chain

    Thanks to:
    DSSP: W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637
    (secondary structure determination)
    dssp2pdb (James Stroud, 2002) - converting DSSP secondary structure
    assignments to pdb format
    """
    outputPDB = []
    #generate pdb for each chain; add TER at the end of each chain
    for chain in splitPDB:
        outputPDB = outputPDB + generatePDB(chain)
        if chain[0].groupType == 'std':
            outputPDB = outputPDB + ['TER\n']

    #Generate secondary structure
    #write out temporary pdb
    writeToPDB(outputPDB, 'tmp.pdb')
    os.system('dssp tmp.pdb tmp.dssp')
    os.system('dssp2pdb tmp.dssp > tmp.pdbheader')
    pdbheader = open('tmp.pdbheader', 'r').readlines()
    os.system('rm -f tmp.pdb tmp.dssp tmp.pdbheader')

    return pdbheader + outputPDB

def writeToPDB(outputLines, filename):
    outputFile = open(filename, 'w')
    for line in outputLines:
        outputFile.write(line)

def writeSinglePDB(loadedPDB, filename):
    #Write a list of PDB lines to a file (from generatePDB)
    #split by chains
    parsedPDB = parseByChain(loadedPDB)
    #generate and write a pdb model, including TERs and secondary structure
    outputPDB = generate_one_model(parsedPDB)
    writeToPDB(outputPDB, filename)

def writeMultipleModels(loadedPDBlist, filename):
    """
    For writing multiple models of the same protein to a pdb file.

    1) Each model must be bracketed by 'MODEL' and 'ENDMDL' keywords so
    that viewer programs do not bond distinct models together
    2) Each chain must have a unique chainID so that they can be
    distinguished by viewer programs
    3) secondary structure and TERs (see generate_one_model)
    """
    outputPDB = []
    chainList = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
                 'P','Q','R','S','T','U','V','W','X','Y','Z']
    chainNum = 0
    for model in loadedPDBlist:
        new_model = []
        for res in model:
            new_model.append(res.copy())
        #add model keyword
        outputPDB = outputPDB + ['MODEL\n']
        parsedModel = parseByChain(new_model)
        #re-chain the model (so each chain in the pdb will be unique)
        for chain in parsedModel:
            if chain[0].groupType == 'std':
                resetChain(chain, chainList[chainNum])
                chainNum = chainNum + 1
        outputPDB = outputPDB + generate_one_model(parsedModel)
        outputPDB = outputPDB + ['ENDMDL\n']
    writeToPDB(outputPDB, filename)


#*******************************************************************************
#Analysis functions for loadedPDBs and PDBres residues

def center_at_origin(pdb):
    newPDB = []
    #center = array(calcCAcentroid(pdb)) LIZA
    center= numarray.array(calcCAcentroid(pdb))
    print "$$$$$$$$$$$$$$$$$$$$$$$$calcCAcentroid", calcCAcentroid(pdb)
    print "$$$$$$$$$$$$$$$$$$$$$$$$center", center

    for res in pdb:
        newRes = res.copy()
        #print newRes.printid
        for atom in newRes.atomlist:
            #print atom
            #print newRes[atom]
            print "$$$$$$$$$$$$$newRes[atom in before after", newRes[atom]
            newRes[atom] = array(newRes[atom])
            print "$$$$$$$$$$$$$newRes[atom in center after", newRes[atom]
            newRes[atom] = newRes[atom] - center
            #print newRes[atom][0]
            #print ''
        newPDB.append(newRes)
    return newPDB

def translatePDB(pdb, vector):
    newPDB = []
    for res in pdb:
        newRes = res.copy()
        for atom in newRes.atomlist:
            #print atom
            #print newRes[atom]
            newRes[atom] = numarray.array(newRes[atom])
            newRes[atom] = newRes[atom] + vector
            #print newRes[atom][0]
            #print ''
        newPDB.append(newRes)
    return newPDB

def getCentroid(atomlist):
    coordlist=[]
    for item in atomlist:
        coordlist.append(item.coords)
    return coordlib.vavg(coordlist)

def addIndices(loadedPDB,start=0):
    """
    Adds numerical indices to residues in a loadedPDB to facilitate
    later data processing.
    """
    Nres = len(loadedPDB)
    for i in range(0, Nres):
        setattr(loadedPDB[i], 'index', i+start)

def popLigands(loadedPDB):
    """
    removes ligands from a loadedPDB object and returns them.
    Input:   a loadedPDB object (output of loadProt)
    Output:  ligands are removed from loadedPDB
             function returns the ligands

    Intended usage:  Some protein structure calculations use ligands (e.g.
    contacts), while others must ignore them (protein torsion angles).
    For the latter kind of calculation, this allows you to remove the ligands
    during the calculation but have them handy to add back in at the end to
    preserve the pdb.
    """
    #Store protein and ligand atoms
    protPDB = []
    ligandPDB = []
    for res in loadedPDB:
        if res.groupType != 'het':
            protPDB.append(res)
        else:
            ligandPDB.append(res)
            #sys.stderr.write('removing hetero group ' + res.printid + '\n')
    #reset loadedPDB (change in place) to protPDB
    loadedPDB[:] = protPDB
    #return the ligands
    return ligandPDB

def getLigands(loadedPDB):
    #return all ligand residues from loadedPDB
    ligands = []
    for res in loadedPDB:
        if res.groupType == 'het':
            ligands.append(res)
    if len(ligands) == 0:
        sys.stderr.write('warning: found no ligands in this pdb\n')
    return ligands

def extractChains(loadedPDB, chains):
    #return any chain in chains (more than one allowed) from loadedPDB
    chainPDB = []
    for res in loadedPDB:
        if res.chain in chains:
            chainPDB.append(res)
    if len(chainPDB) == 0:
        print 'error:  cannot find chain', chains, 'in this PDB'
        print 'Check your chains and try again'
        sys.exit(1)
    return chainPDB

def parseByChain(loadedPDB):
    """
    Parse a loadedPDB into chains.
    Input:  loadedPDB (output of loadProt())
    Output:  a nested list of chains corresponding to loadedPDB

    Just a standard parsing routine.
    """
    #ligands = popLigands(loadedPDB)
    chainList = [[loadedPDB[0]]]
    Nres = len(loadedPDB)
    for i in range(1, Nres):
        if loadedPDB[i].chain == loadedPDB[i-1].chain:
            chainList[-1].append(loadedPDB[i])
        else:
            chainList.append([loadedPDB[i]])
    #if ligands != []:
    #    chainList.append(ligands)
    return chainList

def resetChain(loadedPDB, newChain):
    #Resets the chain of all residues of a loadedPDB to a new chain
    #Don't do anything if the chainID is already correct
    if loadedPDB[0].chain == newChain:
        return
    for res in loadedPDB:
        res.chain = newChain

def combine_identical_chains(chainParsedPDB):
    """
    Input:  chainParsedPDB (output of parseByChain(loadedPDB))
    Method:  using Nres of each chain, detect identical chain.  For ex., if I
    had four chains with following Nres:

    [141,146,141,146]

    I would split the list as follows:

    [ [141,141], [146,146] ]

    This allows me to compare identical monomers
    """
    Nchains = len(chainParsedPDB)
    #make NresList
    NresList = []
    for chain in chainParsedPDB:
        NresList.append(len(chain))
    #Lists to keep track of matches
    matchedList = [0] * Nchains
    combinedList = []
    for i in range(Nchains):
        #skip chains already matched.
        if matchedList[i] == 1:
            continue
        imatchList = [chainParsedPDB[i]]
        for j in range(i+1,Nchains):
            if NresList[i] == NresList[j]:
                imatchList.append(chainParsedPDB[j])
                matchedList[j] = 1
        combinedList.append(imatchList)
        matchedList[i] = 1
    return combinedList

def reset_chain(loadedPDB, chain):
    for res in loadedPDB:
        res.chain = chain
        res.printid = res.restype + ' ' + res.chain + ' ' + res.resid

def extractResidues(loadedPDBchain, resnum1, resnum2):
    """
    inputs:
    loadedPDBchain - a single chain of a pdb
    resnum1, resnum2 - the range to extract

    No protection against multiple chains
    """
    residueList = []
    for res in loadedPDBchain:
        if res.resnum >= resnum1:
            residueList.append(res)
            if res.resnum == resnum2:
                return residueList
    print 'residues', str(resnum1) + '-' + str(resnum2), 'not found'
    sys.exit(1)

def getSequence(loadedPDB):
    seqStr= ''
    for res in loadedPDB:
        if res.groupType == 'het':
            continue
        seqStr = seqStr + res.resaa
    return seqStr

def addasa(loadedPDB, pdbfile):
    """
    Call naccess on a pdb and then add the information into each residue of
    the loadedPDB
    """
    pdbLigands = popLigands(loadedPDB)
    pdbid = string.split(pdbfile, '.')[0]
    os.system('naccess ' + pdbfile)
    rsafname = pdbid + '.rsa'
    rsafile = open(rsafname, 'r')
    i=0
    for line in rsafile:
        if line[0:3] == 'RES':
            res = loadedPDB[i]
            setattr(loadedPDB[i], 'asatot', float(string.strip(line[22:28])))
            setattr(loadedPDB[i], 'asasc', float(string.strip(line[35:41])))
            setattr(loadedPDB[i], 'asabb', float(string.strip(line[48:54])))
            i=i+1
    loadedPDB[:] = loadedPDB + pdbLigands
    os.system('rm -f ' + rsafname)
    os.system('rm -f ' + pdbid + '.asa')

def getCAlist(loadedPDB):
    CAlist = []
    for res in loadedPDB:
        try:
            CAlist.append(res.CA)
        except AttributeError:
            continue
    return CAlist

def calcCAcentroid(loadedPDB):
    #returns centroid of the CA atoms from a loadedPDB
    CAlist = getCAlist(loadedPDB)
    return coordlib.vavg(CAlist)

def get_index(loadedPDB, resnum):
    #returns the index of the loadedPDB residue with residue number resnum.
    Nres = len(loadedPDB)
    for i in range(Nres):
        if loadedPDB[i].resnum >= resnum:
            return i
    if loadedPDB[Nres-1].resnum <= resnum:
        return Nres-1

def transformResCoords(PDBres, tran, rot):
    newRes = PDBres.copy()
    print "$$$$$$$$$rot", rot
    for atom in newRes.atomlist:
    	print "$$$$$$$$$$$$$newRes[atom] before:", newRes[atom]
        newRes[atom] = numarray.array(newRes[atom])
        print "$$$$$$$$$$$$$newRes[atom] after:", newRes[atom]
        newRes[atom] = numarray.matrixmultiply(newRes[atom], rot) + tran
        print "$$$$$$$$$$$$$newRes[atom] after matrixmult:", newRes[atom]
    return newRes

def transformCoords(loadedPDB, tran, rot):
    """
    Input:  loadedPDB, a translation vector, and a rotation matrix
    Output: a new loadedPDB with the desired transformations
    """
    newPDB = []
    for res in loadedPDB:
        newRes = transformResCoords(res, tran, rot)
        newPDB.append(newRes)
    return newPDB

def axisangle_to_rot(axis, angle):
    """
    Converts a rotation axis and angle to a rotation matrix where
    a clockwise rotation is a positive value of 'angle'
    from Martin Baker
    http://www.euclideanspace.com/maths/geometry/rotations/conversions/
    angleToMatrix/index.htm
    """
    c,s = cos(angle), sin(angle)
    t = 1-c
    x,y,z = axis[0], axis[1], axis[2]

    rot = [[t*x*x + c, t*x*y - z*s, t*x*z + y*s],
     [t*x*y + z*s, t*y*y + c, t*y*z - x*s],
     [t*x*z - y*s, t*y*z + x*s, t*z*z + c]]
    return rot

def condenseReslist(resList):
    if resList == []:  return []
    rasmol_reslist = ''
    resRangeList = [[resList[0]]]
    Nres = len(resList)
    for i in range(1, Nres):
        if int(resList[i].resid) - int(resList[i-1].resid) == 1:
            continue
        resRangeList[-1].append(resList[i-1])
        resRangeList.append([resList[i]])
    resRangeList[-1].append(resList[-1])
    Nrange = len(resRangeList)

    condensed_reslist = []
    for i in range(Nrange):
        res = resRangeList[i]
        if res[0].resid == res[1].resid:
            condensed_reslist.append(res[1].resid + res[1].chain)
        else:
            condensed_reslist.append(res[0].resid + '-' + res[1].resid +
                                     res[1].chain)
    return condensed_reslist

def printResList(resList):
    if resList is None:
        return 'None'
    outstr = ''
    condensedResList = condenseReslist(resList)
    Nres = len(condensedResList)
    for i in range(Nres):
        outstr = outstr + condensedResList[i]
        if i < Nres - 1:
            outstr = outstr + ','
    return outstr
