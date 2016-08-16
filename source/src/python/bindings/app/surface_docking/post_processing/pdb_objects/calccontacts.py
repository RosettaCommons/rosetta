#!/usr/bin/env python

from loadPDB import *
import coordlib
import string
import math
import copy

class atom_contact:
    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2
        #hbonds only relevant to full atom mode
        if contactType == 'full':
            self.isHbond = (self.atom1[0] in ['N', 'O'] and self.atom2[0]
                            in ['N', 'O'])

    def numSwap(self):
        #swap the identities of atom1 and atom2
        tmp_atom = self.atom1
        self.atom1 = self.atom2
        self.atom2 = tmp_atom

class resContact:
    def __init__(self, partnerA_res, partnerB_res, atomContactList):
        self.res1 = partnerA_res
        self.res2 = partnerB_res
        self.index1 = partnerA_res.index
        self.index2 = partnerB_res.index
        self.atomContactList = atomContactList
        #only relevant to full atom mode
        if contactType == 'full':
            self.Ncontacts = len(atomContactList)
            self.Nhbonds = 0
            for atom_contact in atomContactList:
                if atom_contact.isHbond:
                    self.Nhbonds = self.Nhbonds + 1

    def copy(self):
        tmp_copy = copy.copy(self)
        tmp_copy.atomContactList = self.atomContactList[:]
        Natomcontacts = len(self.atomContactList)
        for i in range(Natomcontacts):
            tmp_copy.atomContactList[i] = atom_contact \
            (self.atomContactList[i].atom1, self.atomContactList[i].atom2)
        return tmp_copy

    def reduceContact(self, res):
        """
        This changes the form of resContact from a contact-map type format
        to a contact-profile type format.  A contact map is a list of
        contacts, whereas a contact profile is a list of residues with
        their contacts.

        In a contact profile, each residue only needs to know which other
        residues it contacts, not the complete resContact format.

        When creating a contact profile, each contact in the map is added
        to both residues involved.  This means changing the format from:

        i-j
        i-k
        i-l
        j-k
        j-l
        j-m

        (list/map format)

        to

        i-j,k,l
        j-i,k,l,m

        (profile format)

        For example, for residue i, I want to keep only the info for j,k, and
        l since I already have the info for residue i

        for residue j, I want to keep only the information for k,l, and m.  For
        the i-j contact it is more complicated because I need to invert the
        order to j-i before removing residue j.

        Ex. for i-j contact

        1) make one copy
        2) delete the info for residue i ('res1', 'index1')
        3) add contact to residue i

        4) make another copy
        5a) invert identity of res1/res2, index1/index2
         b) invert each of the atom contacts
        6) delete res1, index1 which now match j.

        Input parameter res determines which option is picked
        If 'res1', carry out steps 1-3
        If 'res2', carry out steps 4-6.
        """
        if res == 'res2':
            #swap residue indices
            self.res2 = self.res1
            self.index2 = self.index1
            #swap atom contacts
            for atom_contact in self.atomContactList:
                atom_contact.numSwap()
        delattr(self, 'res1')
        delattr(self, 'index1')

class shellAtom:
    def __init__(self, printid, index, atomname):
        """
        Class attributes:

        printid:  printid of the atom's residue
        index:    index of the atoms residue
        atomname:  the atom name
        count:  number of times that shellAtom occurs (for multiple monomers)
        atomid:  a tag for easy comparison of shell atoms
        """
        self.printid = printid
        self.index = index
        self.atomname = atomname
        self.count = 1
        self.atomid = str(self.index) + self.atomname

def addMax_radius(loadedPDBres):
    #calculate maximum radius
    if loadedPDBres.max_radius != None:
        return
    max_radius_sq = 0
    for atom in loadedPDBres.coordlist:
        atomdist_sq = coordlib.dist_sq(atom, loadedPDBres.centroid)
        max_radius_sq = max(atomdist_sq, max_radius_sq)
    setattr(loadedPDBres, 'max_radius', math.sqrt(max_radius_sq))

def makeContactPDB(loadedPDB):
    """
    Add certain attributes to loadedPDB to make it a contactPDB:
    centroid:  centroid of the residue or ligand
    max_radius:  farthest atom from the centroid

    To be computed later:
    contactList:  a list of all contacts made by that residue
    Ncontacts
    Nhbonds
    """
    for res in loadedPDB:
        #initialize centroid only in full atom mode
        if getAtoms == 'centroid':
            res.addCentroid()
            addMax_radius(res)
        #initialize sccentroid only in sccentroid mode
        if getAtoms == 'sccentroid':
            res.addscCentroid()
        #Ncontacts and Nhbonds should only be set in full atom mode
        if contactType == 'full':
            setattr(res, 'Ncontacts', 0)
            setattr(res, 'Nhbonds', 0)
        setattr(res, 'contactList', [])
        #NresContacts (res-res pairs) is relevant to all modes
        setattr(res, 'NresContacts', 0)

def parse_and_pair(contactPDB, parserule):
    """
    This is the major function that allows one find-contacts script to
    calculate contacts in a variety of different ways.
    Using the parserule, the appropriate segment(s) of the pdb are extracted.

    Two types of contact calculations are enabled by this routine, and the
    type of calculation determines the inputs and actions:
    internal contacts   (no '-' in the parserule)
    interface/binding site contacts  ('-' in the parserule determining which
    subsets of the protein to use)

    For internal calculations:
    The following types of input are accepted:
    'all' - all residues in the protein, including all chains and any ligand
    groups in the pdb
    'A' - chain A only, all internal contacts
    'AB' - chains A and B, all internal + A-B interface

    action:  this subroutine extracts the appropriate portion of the pdb and
    then passes it to pair_close_internal(contactPDB), and returns the result

    For interface calculations:
    'A-B', 'AB-C', 'AB-CD', etc:  docking partners.  Can be multiple chains
    per partner
    'A-ligand', 'AB-ligand', 'ligand-A' - subset of protein vs. all ligands
    'all-ligand' - all protein vs. all ligands
    action:  extracts the units separated by '-' and passes them through
    pair_close_interface, and returns the result.


    Input:  contactPDB (see class contactProtein)
    Output:  [partnerA, partnerB] where partnerA and partnerB are subsets of
    the contactPDB such that partnerA-partnerB contacts are calculated but not
    internal contacts within either partner A or partner B.  A parsing rule is
    used to split loadePDB into partners A and B.
    Currently, 3 parsing rules are implemented:

    This routine is flexible so that new parsing rules can be
    added independently of any previous rules
    """
    #interface contacts mode
    if '-' in parserule:
        splitRule = string.split(parserule, '-')
        #protein-ligand contacts
        if 'ligand' in splitRule:
            # 'ligand-A' form
            if splitRule[0] == 'ligand':
                proteinChain = splitRule[1]
            # 'A-ligand' form
            else:
                proteinChain = splitRule[0]
            # 'all-ligand'
            # Pop ligands from contactPDB and assign to partner B
            partnerB = popLigands(contactPDB)
            if proteinChain == 'all':
                partnerA = contactPDB
            else:
                partnerA = extractChains(contactPDB, proteinChain)
        # two docking/interface partners (e.g. 'AB-C')
        else:
            partnerA = extractChains(contactPDB, splitRule[0])
            partnerB = extractChains(contactPDB, splitRule[1])
        #see pair_close_interface subroutine
        return pair_close_interface(partnerA, partnerB)
    #internal contacts
    else:
        #whole PDB
        if parserule == 'all':
            subset = contactPDB
        #some subset (specified chains)
        else:
            subset = extractChains(contactPDB, parserule)
        #see pair_close_internal subroutine below
        return pair_close_internal(subset)

def pair_close_internal(contactPDBsubset):
    """
    Input:  subset of contact PDB

    Output:  all pairs of residue pairs i, j where a contact is possible
    according to the reference (max radius to centroid) distances for
    each group plus 5.0A
    """
    Nres = len(contactPDBsubset)
    close_pair_list = []
    for i in range(0, Nres):
        for j in range(i+1, Nres):
            #refChain filter
            if refChain != None:
                if not(contactPDBsubset[i].chain in refChain or
                       contactPDBsubset[j].chain in refChain):
                    continue

            #max_radius from centroid
            #filter = max_radius1 + max_radius2 + 5.0A
            #if the centroids are farther apart than this, there can't be any
            #contacts between residues i and j
            thispair = 0
            #pair without filtering (if global flag smartfilter is not turned
            #on)
            if not(smartfilter):
                thispair = 1
            #pair with filtering on centroid - centroid distance (applicable
            #only to full atom mode)
            else:

                filter_sq = (contactPDBsubset[i].max_radius +
                             contactPDBsubset[j].max_radius + contactFilter +
                             1.0)**2
                cendist_sq = coordlib.dist_sq(contactPDBsubset[i].centroid,
                                              contactPDBsubset[j].centroid)
                if cendist_sq <= filter_sq:
                    thispair = 1

            if thispair:
                close_pair_list.append([contactPDBsubset[i],
                                        contactPDBsubset[j]])
    return close_pair_list

def pair_close_interface(partnerA, partnerB):
    """
    Input:  partner A and partner B of a contactPDB
    In internal case, just [partnerA] is given to save memory

    Output:  all pairs of residues between partner A and
    partner B for which a contact is possible according to the reference
    distances for each group
    """
    close_pair_list = []
    NresA = len(partnerA)
    NresB = len(partnerB)
    for i in range(0, NresA):
        for j in range(0, NresB):
            #max_radius from centroid
            #filter = max_radius1 + max_radius2 + 5.0A
            #if the centroids are farther apart than this, there can't be any
            #contacts between residues i and j
            thispair = 0
            #pair without filtering (if global flag smartfilter is not turned
            #on)
            if not(smartfilter):
                thispair = 1
            #pair with filtering on centroid - centroid distance (applicable
            #only to full atom mode)
            else:
                filter_sq = (partnerA[i].max_radius + partnerB[j].max_radius +
                             contactFilter + 1.0)**2
                cendist_sq = coordlib.dist_sq(partnerA[i].centroid,
                                              partnerB[j].centroid)
                if cendist_sq <= filter_sq:
                    thispair = 1
            if thispair:
                close_pair_list.append([partnerA[i], partnerB[j]])
    return close_pair_list

def filterContacts(close_residue_list):
    """
    Filter a close_residue_list (output of pair_close_refAtoms) for residues
    with at least one atom-atom contact of < 4.0A+
    (the filter can be set to another distance in class contactProtein)
    """

    resContactList = []
    for close_pair in close_residue_list:
        makesContact = 0
        partnerA_res = close_pair[0]
        partnerB_res = close_pair[1]
        #bond filter - do not count covalent bonds between residues (< 2.0A)
        #as contacts
        bondfilter_sq = 4.0
        contactfilter_sq = contactFilter**2
        atomContactList=[]
        #full atom mode; iterate over all atoms
        if contactAtoms == 'all':
            idA = int(partnerA_res.resid)
            idB = int(partnerB_res.resid)
            A_atoms = partnerA_res.atomlist
            B_atoms = partnerB_res.atomlist
        #designed if you only want contacts between CA-CA, centroid-centroid,
        #etc.
        #Just count contacts for those atom pairs
        else:
            A_atoms = B_atoms = [contactAtoms]
        #all atoms from partner A res against all atoms of partner B res
        for A_atom in A_atoms:
            for B_atom in B_atoms:
                atomdist_sq = coordlib.dist_sq(partnerA_res[A_atom],
                                               partnerB_res[B_atom])
                if atomdist_sq >= bondfilter_sq and atomdist_sq > 16.0 and atomdist_sq <= 64.0: # kepler
                #if atomdist_sq >= bondfilter_sq and atomdist_sq <= contactfilter_sq: #Mike
                    #see atom_contact class above
                    atomContactList.append(atom_contact(A_atom, B_atom))
                #if idA > 23 and idA < 33 and idB > 23 and idB < 33 and idB - idA > 5:
                #    if atomdist_sq >= 16.0 and atomdist_sq <= 36.0:
                #        if A_atom[:1] == 'C' or A_atom[:1] == 'N':
                #            if B_atom[:1] == 'C' or B_atom[:1] == 'N':
                #                print str(idA)+A_atom+str(idB)+B_atom

        if atomContactList != []:
            #see resContactList class above
            resContactList.append(resContact(partnerA_res, partnerB_res,
                                             atomContactList))

    return resContactList

def makeContactProfile(contactMap, loadedPDB):
    """
    Takes the contactMap of a contactProtein structure and adds the contacts
    into the original loadedPDB
    In makeContactPDB(loadedPDB), an empty contactList was added to each
    residue in the loaded PDB.  This list will now be filled up with all the
    contacts that each residue makes.

    Also fill in Ncontacts and Nhbonds for each residue, which were initialized
    to zero.
    """
    #iterate through the contact map
    #The loadedPDB indices of each contact have been recorded and will now
    #be used to send each contact to the contact lists of the appropriate
    #two residues
    for contact in contactMap:
        res1 = loadedPDB[contact.index1]
        res2 = loadedPDB[contact.index2]
        #Make two copies of the contact, one for res1 and one for res2.
        #Reset the resContact and add it to the appropriate residue.
        contact_copy1 = contact.copy()
        contact_copy2 = contact.copy()
        contact_copy1.reduceContact('res1')
        contact_copy2.reduceContact('res2')
        res1.contactList.append(contact_copy1)
        res2.contactList.append(contact_copy2)
        res1.NresContacts = res1.NresContacts + 1
        res2.NresContacts = res2.NresContacts + 1
        #only relevant to full atom mode!
        if contactType == 'full':
            res1.Ncontacts = res1.Ncontacts + contact.Ncontacts
            res1.Nhbonds = res1.Nhbonds + contact.Nhbonds
            res2.Ncontacts = res2.Ncontacts + contact.Ncontacts
            res2.Nhbonds = res2.Nhbonds + contact.Nhbonds

def makeShellProfile(contactProfile):
    """
    This converts a contact profile (list of residue-residue contacts for each
    residue) to a shell profile (list of the atomic shell for each residue).

    Contact profile format is useful for listing residues that contact a
    particular residue, but it is tedious for analyzing the atomic shell of a
    residue since one must search the atoms of each residue in the list.  Here
    is an example:

    contact profile format:

    contacts for residue LYS A 4:
    GLU A 3  N,CA,OE1
    PRO A 5  CA,CB,CG

    contact shell format:

    shell for residue LYS A 4:
    atom 1:  GLU A 3 N
    atom 2:  GLU A 3 CA
    etc.

    This allows direct iteration over the atomic shell without having to iterate
    through particular residues first.

    Practically, for each atom in the shell, I want to keep track of printid
    (residue label), index (list identifier used to compare different pdbs),
    atom name, counts (used if multiple monomers are added), and an
    I will keep an atom label that combines the index and the name for easy
    referencing.  These attributes will be recorded in the shellatom class
    at the top.
    """
    for res in contactProfile:
        setattr(res, 'atomShell', [])
        for resContact in res.contactList:
            atomList = []
            for atom_contact in resContact.atomContactList:
                #record printid, index, and atom name
                if not(atom_contact.atom2 in atomList):
                    #use res.name (the name I want to use for comparisons,
                    #not the printid taken directly from the pdb)
                    atomLabel = resContact.res2.name + atom_contact.atom2
                    #res.atomShell.append(shellAtom(resContact.res2.printid,
                    #resContact.index2, atom_contact.atom2))
                    res.atomShell.append(atomLabel)
                    atomList.append(atom_contact.atom2)
        delattr(res, 'contactList')

smartfilter = 1
contactAtoms = 'all'
contactFilter = 8.0 #kepler
#contactFilter = 4.0 # Mike
contactType = 'full'
getAtoms = 'centroid'
ligandFlag = 1

#This variable can be set externally.  Basically, it limits contact calculation
#to those pairs where at least one residue is in one chain in the string
#refChain.  Is set to None (unrestricted) by default
#(only applies to internal contacts calculations)
refChain = None

def initialize_global_flags(mode):
    """
    global flags dependent on mode:

    smartfilter:  filter on centroid-centroid distance in full atom mode,
    but not in CA-CA mode or sidechain-sidechain centroid mode

    contactAtoms:  defines atom list to iterate over when calculating contacts.
    Allows full atom mode and single atom modes to use the same contact-
    counting function

    contactType:  control creation and printing of certain contact attributes
    (Ncontacts, Nhbonds) only relevant to full atom mode

    contactFilter:  the contact distance cutoff is dependent on the mode.
    Currently, it is 4.0A for full atom mode, 7.0A for CA-CA mode, and 6.0A for
    sidechain-sidechain centroid mode.

    getAtoms:  global flag to determine if (whole residue) centroid and
    sidechain centroid need to be initialized.  To save computing time.
    """
    global smartfilter, contactType, contactAtoms, contactFilter, getAtoms
    if mode == 'fullatom':
        return
    if mode == 'shell':
        contactFilter = 6.0
    if mode in ['calpha', 'sccentroid', 'disulf']:
        smartfilter = 0
        contactType = 'brief'
    if mode == 'calpha':
        contactAtoms = 'CA'
        contactFilter = 7.0
        return
    if mode == 'sccentroid':
        contactAtoms = 'sccentroid'
        contactFilter = 6.0
        getAtoms = 'sccentroid'
        return
    if mode == 'disulf':
        contactAtoms = 'SG'
        contactFilter = 3.0

def addContacts(startProtein, mode, parserule):
    """
    Add contacts to a Protein object from loadPDB.py
    This is somewhat complicated because I need to convert the Protein into a
    loadedPDB to calculate the contacts and then convert it back again.
    However, this works efficiently because the getAllRes() function of a
    Protein object maintains the original order of residues so that when
    the protein is recreated using parse(), the original structure is
    maintained.
    """
    #Create a contact object and calculate contacts
    if ligandFlag:
        newLoadedPDB = startProtein.getAllRes()
    else:
        newLoadedPDB = startProtein.getAllProteinRes()
        print 'ligand flag off; not counting ligands in contacts'
    newContactProtein = contactProtein(mode)
    newContactProtein.contactPDB = newLoadedPDB
    newContactProtein.find_contacts(parserule)
    newContactProtein.initContactProfile()

    #Re-create the protein structure
    newProtein = Protein()
    newProtein.parse(newContactProtein.contactPDB)
    startProtein = newProtein

def makeContactMap(pdb1, pdb2):
    """
    Calculate contacts between any two pdbs
    """
    initialize_global_flags('fullatom')
    makeContactPDB(pdb1)
    makeContactPDB(pdb2)
    closeResiduePairs = pair_close_interface(pdb1, pdb2)
    return filterContacts(closeResiduePairs)

class contactProtein:
    """
    The contactProtein class contains information about the protein's contacts
    in two forms:
    list of contacts
    contact profile (contacts by residue)
    """
    def __init__(self, mode):
        self.mode = mode
        initialize_global_flags(mode)
        self.contactPDB = None

    def initializePDB(self, PDBfile):
        self.contactPDB = loadProt(PDBfile)

    def find_contacts(self, parserule):
        #mode currently allows full atom, CA-CA, and sidechain centroid-
        #centroid
        #set global control flags to control behavior of the functions with
        #respect to mode
        if self.contactPDB == None:
            print 'You are trying to calculate contacts without a loadedPDB'
            print 'Call initializePDB(PDBfile) to load the PDB from a PDB'
            print 'file or load the PDB externally and add it to the \
            contactProtein'
            sys.exit(1)
        if contactType == 'brief':
            popLigands(self.contactPDB)
        makeContactPDB(self.contactPDB)
        addIndices(self.contactPDB)
        closeResiduePairs = parse_and_pair(self.contactPDB, parserule)
        self.contactMap = filterContacts(closeResiduePairs)

    def initContactProfile(self):
        #make a contact profile only if needed (contact map always needed)
        makeContactProfile(self.contactMap, self.contactPDB)
        if self.mode == 'shell':
            self.convert_to_shell()

    def convert_to_shell(self):
        """
        convert from contact profile format (list of residue contacts for each
        residue) to shell format (list of atoms in the atomic shell of a given
        residue)
        """
        makeShellProfile(self.contactPDB)

    def outputMap(self, dest, PDBfile):
        """
        Output method to either print contacts or write them to a file
        dest = 'print' outputs the contact map to the screen
        dest = 'file' outputs the contact map to a file
        """
        if dest == 'file':
            pdbid = string.split(PDBfile, '.')[0]
            outfname = pdbid + '.contactmap'
            output = open(outfname, 'w')
        #header = 'restype1 chain1 resnum1 restype2 chain2 resnum2'
        #if contactType == 'full':
        #    header = header + ' Ncontacts Nhbonds'
        #if dest == 'print':
        #    print header
        #else:
        #    output.write(header + '\n')
        for resContact in self.contactMap:
            printList = [resContact.res1.printid, resContact.res2.printid]
        #    if contactType == 'full':
        #        printList = printList + [resContact.Ncontacts,
       #                                  resContact.Nhbonds]
            outstr = ''
            for item in printList:
                outstr = outstr + str(item) + ' '
            if dest == 'print':
                print outstr
            else:
                output.write(outstr + '\n')

    def outputProfile(self, dest, PDBfile):
        """
        Output method to print contact profile or write them to a file.
        Number of contacts and hbonds for each residue
        equivalent to find_interface.pl in interface mode
        dest = 'print' sends the output to the screen
        dest = 'file sends the output to a file (pdb.contactprofile)
        """
        print 'wtf!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        if dest == 'file':
            pdbid = string.split(PDBfile, '.')[0]
            outfname = pdbid + '.contactprofile'
            output = open(outfname, 'w')
        header = 'restype chain resnum NresContacts'
        if contactType == 'full':
            header = header + ' Ncontacts Nhbonds'
        if dest == 'print':
            print header
        else:
            output.write(header + '\n')
        for res in self.contactPDB:
            if res.NresContacts > 0:
                printList = [res.printid, res.NresContacts]
                if contactType == 'full':
                    printList = printList + [res.Ncontacts, res.Nhbonds]
                outstr = ''
                for item in printList:
                    outstr = outstr + str(item) + ' '
                if dest == 'print':
                    print outstr
                else:
                    output.write(outstr + '\n')

