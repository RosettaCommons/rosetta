#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

'''A simple PDB reader/writer for basic manipulations.

Each line is preserved as much as possible, so that i/o is lossless.
I deliberately modeled the PDB file as a list of lines rather than some more
complex object so that it would be easy to slice up using for comprehensions:
that's why the functions are not member functions of some "PdbFile" object.
'''

import math

class PdbRecord:
    '''A generic line from a PDB file'''
    # line          80+ char string without trailing newline

    def __init__(self, line):
        d = self.__dict__ # direct assignment may cause endless loops
        d['line'] = line.rstrip("\n").ljust(80)

    def __setattr__(self, name, value):
        '''String representation is automatically updated when a value is set'''
        if name == 'line': raise ValueError("cannot change value of 'line' directly")
        self.__dict__[name] = value

    def __str__(self):
        return self.line

class Atom (PdbRecord):
    '''An ATOM or HETATM line from a PDB file'''
    # het           boolean, is this a HETATM?
    # name          4-character atom name e.g. " CA "
    # resname       3-character residue type name e.g. "ALA"
    # chain         single-character chain ID
    # resnum        residue sequence number e.g. 42
    # inscode       residue insertion code e.g. " "
    # resseq        residue number + insertion code, as a 5-character string
    # x, y, z       floating-point coordinates of the atom/hetatm
    # segid         PDB "segment identifier", mostly unused (73-76)

    def __init__(self, line):
        PdbRecord.__init__(self, line)
        d = self.__dict__ # direct assignment may cause endless loops
        d['het']        = self.line.startswith("HETATM")
        d['name']       = self.line[12:16]
        d['resname']    = self.line[17:20]
        d['chain']      = self.line[21:22]
        d['resnum']     = int(self.line[22:26])
        d['inscode']    = self.line[26:27]
        d['resseq']     = self.line[22:27]
        d['x']          = float(self.line[30:38])
        d['y']          = float(self.line[38:46])
        d['z']          = float(self.line[46:54])
        d['segid']      = self.line[72:76]

    def __setattr__(self, name, value):
        '''String representation is automatically updated when a value is set'''
        PdbRecord.__setattr__(self, name, value)
        l = self.line # direct assignment may cause endless loops
        if name == 'het':
            if value:           l = "HETATM" + self.line[6:]
            else:               l = "ATOM  " + self.line[6:]
        elif name == 'name':    l = "%s%4.4s%s" % (self.line[:12], value, self.line[16:])
        elif name == 'resname': l = "%s%3.3s%s" % (self.line[:17], value, self.line[20:])
        elif name == 'chain':   l = "%s%1.1s%s" % (self.line[:21], value, self.line[22:])
        elif name == 'resnum':  l = "%s%4.4s%s" % (self.line[:22], value, self.line[26:])
        elif name == 'inscode': l = "%s%1.1s%s" % (self.line[:26], value, self.line[27:])
        elif name == 'resseq':  l = "%s%5.5s%s" % (self.line[:22], value, self.line[27:])
        elif name == 'x':       l = "%s%8.3f%s" % (self.line[:30], value, self.line[38:])
        elif name == 'y':       l = "%s%8.3f%s" % (self.line[:38], value, self.line[46:])
        elif name == 'z':       l = "%s%8.3f%s" % (self.line[:46], value, self.line[54:])
        elif name == 'segid':   l = "%s%4.4s%s" % (self.line[:72], value, self.line[76:])
        self.__dict__['line'] = l

    def is_H(self):
        '''Returns true iff this atoms looks like a hydrogen'''
        c1 = self.name[0]
        c2 = self.name[1]
        return ((c1 == 'H' and c2 != 'G') or (c1.isdigit() and c2 == 'H'))

def read_pdb_file(f):
    '''Returns a list of PdbRecord objects.

    f is either a file name or a file object.
    '''
    if isinstance(f, str):
        f = open(f, 'rU')
        ret = read_pdb_file(f)
        f.close()
        return ret
    records = []
    for line in f:
        card = line[:6]
        if   card == 'ATOM  ': records.append(Atom(line))
        elif card == 'HETATM': records.append(Atom(line))
        else: records.append(PdbRecord(line))
    return records

def write_pdb_file(f, records):
    '''Writes a list of PdbRecord objects to file.

    f is either a file name or a file object.
    '''
    if isinstance(f, str):
        f = open(f, 'rU')
        ret = write_pdb_file(f, records)
        f.close()
        return ret
    for r in records:
        f.write(r.line)
        f.write("\n")

def get_atoms(records):
    '''Simple function to fish out all ATOM and HETATM records.'''
    return [r for r in records if isinstance(r, Atom)]

def get_het_atoms(records):
    '''Simple function to fish out all HETATM records.'''
    return [r for r in records if isinstance(r, Atom) and r.het]

def centroid(atoms):
    '''Returns an unweighted center-of-mass for the list of Atom objects as a tuple (x,y,z).'''
    x = 0.0; y = 0.0; z = 0.0
    for atom in atoms:
        x += atom.x
        y += atom.y
        z += atom.z
    l = len(atoms)
    return (x/l, y/l, z/l)

def translate(atoms, x, y, z):
    '''Translates the Atoms by the specified amounts.'''
    for atom in atoms:
        atom.x += x
        atom.y += y
        atom.z += z

def rmsd(atoms1, atoms2):
    '''Calculates root-mean-square deviation without superposition'''
    if len(atoms1) != len(atoms2):
        raise ValueError("Atom lists must be same length")
    d = 0.0
    for a1, a2 in zip(atoms1, atoms2):
        d += (a1.x-a2.x)**2 + (a1.y-a2.y)**2 + (a1.z-a2.z)**2
    return math.sqrt( d / len(atoms1) )

def bounding_box(atoms):
    '''Calculates a simple bounding box as (xmin,xmax,ymin,ymax,zmin,zmax).'''
    xs = [a.x for a in atoms]
    ys = [a.y for a in atoms]
    zs = [a.z for a in atoms]
    return (min(xs), max(xs), min(ys), max(ys), min(zs), max(zs))

def split_chains(records):
    '''Returns a list of lists of Atom objects, one for each chain.

    Non-Atom records are not returned but are used to find chain breaks.
    We take TER, END, and change of existing chain ID as signals of chain break.
    Distances between residues are NOT checked.
    We assume that if multiple models exist, they have already been split.
    '''
    chains = []
    chain = []
    idx = 0
    last_chain = None
    for record in records:
        if isinstance(record, Atom):
            if last_chain == None: last_chain = record.chain
            elif last_chain != record.chain:
                idx += 1
                last_chain = record.chain
                if len(chain) > 0:
                    chains.append(chain)
                    chain = []
            chain.append(record)
        else:
            card = record.line[0:6].rstrip()
            if card == 'TER' or card == 'END':
                idx += 1
                last_chain = None # or we get double increment
                if len(chain) > 0:
                    chains.append(chain)
                    chain = []
    if len(chain) > 0:
        chains.append(chain)
    return chains

def uniquify_chains(records):
    '''Creates unique chain names (within each model) without regard to the original names.

    Should operate on all PDB records, not just Atoms.
    We take TER, END, and change of existing chain ID as signals of chain break.
    Distances between residues are NOT checked.

    TODO: This should be re-written to use split_models() and split_chains().
    '''
    chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890abcdefghijklmnopqrstuvwxyz"
    idx = 0
    last_chain = None
    for record in records:
        if isinstance(record, Atom):
            if last_chain == None: last_chain = record.chain
            elif last_chain != record.chain:
                last_chain = record.chain
                idx += 1
            record.chain = chains[idx]
        else:
            card = record.line[0:6].rstrip()
            if card == 'MODEL' or card == 'ENDMDL':
                idx = 0 # reset
                last_chain = None
            elif card == 'TER' or card == 'END':
                idx += 1
                last_chain = None # or we get double increment

'''Dictionary of Rosetta atom names to plausible PDB atom names (element only)'''
rosetta_to_pdb_names = {'CNH2':' C  ', 'COO ':' C  ', 'CH1 ':' C  ',
    'CH2 ':' C  ', 'CH3 ':' C  ', 'aroC':' C  ', 'Ntrp':' N  ', 'Nhis':' N  ',
    'NH2O':' N  ', 'Nlys':' N  ', 'Narg':' N  ', 'Npro':' N  ', 'OH  ':' O  ',
    'ONH2':' O  ', 'OOC ':' O  ', 'S   ':' S  ', 'Nbb ':' N  ', 'CAbb':' C  ',
    'CObb':' C  ', 'OCbb':' O  ', 'Hpol':' H  ', 'Hapo':' H  ', 'Haro':' H  ',
    'HNbb':' H  ', 'Phos':' P  '}

def standardize_atom_names(atoms, move_to_segid=True):
    '''Converts Rosetta atom names to more normal PDB-like names.

    If move_to_segid is True, the old names will overwrite the atom segid field.'''
    for atom in atoms:
        if atom.name in rosetta_to_pdb_names:
            if move_to_segid:
                atom.segid = atom.name
            # Dumb way, doesn't uniquify atom names:
            atom.name = rosetta_to_pdb_names[ atom.name ]

def sort_H_to_end(atoms):
    '''In-place stable sort of hydrogens after all heavy atoms, like Rosetta'''
    atoms.sort(lambda a,b: cmp(a.is_H(), b.is_H()))
