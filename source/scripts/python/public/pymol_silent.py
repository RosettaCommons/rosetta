#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   pymol_silent.py
## @brief  Utility file for silent file loading from, e.g. Pymol
## @author Rocco Moretti (rmoretti@u.washington.edu)

"""pymol_silent.py - a utility program to read Rosetta Silent files externally.

The intention of this script is to allow direct reading of Rosetta formatted "binary" silent files from pymol.
(As a side effect, it can allow a limited extraction of PDBs from binary silent files from the command line.)

Usage:
------

To use, load the file in Pymol (for frequent use, add command to your .pymolrc file:

PyMOL> run rosetta/rosetta_source/src/python/apps/public/pymol_silent.py

Once loaded, you can examine and load files from a binary-formatted silent file:

PyMOL> silent_tags path/to/file/silentfile.out
# List of availible tags printed here
PyMOL> silent_load path/to/file/silentfile.out, [tag1, tag2, ...]

Omiting the tags will load all the ones in the file.

Note that silent_load re-reads the full silent file on each invocation.
For speed reasons, it behooves you to specify as many tags per invocation,
rather than having multiple invocations with a small number of tags.

To assist in efficient loading, the tags are actually python style regexes.
(Set up so that they must match on the entire tag). If you put in a tag as-is,
it should match itself, but more clever expressions are availible for those who
know the sytax. (e.g. "1QYS_0003_00[12]." will pull out all the structures with
the second nstruct number in the 0010 to 0029 range.)

Caveat/Warnings/Limitations:
----------------------------

This script only works with binary formatted silent files. Given the differences in storage format,
it's highly unlikely that it will ever be extended to protein or atom tree diff silent files.

The binary silent file format is highly dependant on the internal representation of datastructures
used by the producing application. (Specifically, how it represents a vector of a class with three
single-precision floats as member variables.) This script assumes the most likely representation,
but if your compiler does funky stuff with padding or layout or floating point representation,
it might not work properly. (But in this case your binary silent file probably won't be transferable
to another machinery portable anyway ...) That said, there is a limited amount of support for
big-endian/little endian autodetection (about the same as in standard Rosetta).

In order to avoid a complete re-write of the database parser (and Rosetta's funky atom reoganization logic),
only a limited number of residue types are fully supported. In particular, arbitrary patches and patch
combinations are not supported. This script should do a reasonable job of showing *something*, although possibly
not with the most asthetically pleasing PyMOL representation.

"""

try:
    from pymol import cmd
except ImportError:
    cmd = None

import sys
import struct
import re

class ResidueDatabase:
    ONE2THREE = {
        'A':'ALA',
        'R':'ARG',
        'N':'ASN',
        'D':'ASP',
        'C':'CYS',
        'Q':'GLN',
        'E':'GLU',
        'G':'GLY',
        'H':'HIS',
        'I':'ILE',
        'L':'LEU',
        'K':'LYS',
        'M':'MET',
        'F':'PHE',
        'P':'PRO',
        'S':'SER',
        'T':'THR',
        'W':'TRP',
        'Y':'TYR',
        'V':'VAL',
    }

    def __init__(self):
        self.atom_names = {}
        self.standard = {}

        # Asterisk means "don't output this atom"
        self.add_standard('ALA', ( 'N', 'CA', 'C', 'O', 'CB', 'H', 'HA', '1HB ', '2HB ', '3HB ',) )
        self.add_standard('ARG', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2', 'H', 'HA', '1HB ', '2HB ', '1HG ', '2HG ', '1HD ', '2HD ', 'HE', '1HH1', '2HH1', '1HH2', '2HH2', ) )
        self.add_standard('ASN', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2', 'H', 'HA', '1HB ', '2HB ', '1HD2', '2HD2', ) )
        self.add_standard('ASP', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2', 'H', 'HA', '1HB ', '2HB ',) )
        self.add_standard('CYD', ( 'N', 'CA', 'C', 'O', 'CB', 'SG', 'H', 'HA', '1HB ', '2HB ',) )
        self.add_standard('CYS', ( 'N', 'CA', 'C', 'O', 'CB', 'SG', 'H', 'HA', '1HB ', '2HB ', 'HG',) )
        self.add_standard('GLN', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2', 'H', 'HA', '1HB ', '2HB ', '1HG ', '2HG ', '1HE2', '2HE2', ) )
        self.add_standard('GLU', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2', 'H', 'HA', '1HB ', '2HB ', '1HG ', '2HG ',) )
        self.add_standard('GLY', ( 'N', 'CA', 'C', 'O', 'H', '1HA ', '2HA ',) )
        self.add_standard('HIS', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2', 'H', 'HA', '1HB ', '2HB ', 'HD2', 'HE1', 'HE2', ) )
        self.add_standard('HIS_D', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2', 'H', 'HA', '1HB ', '2HB ', 'HD1', 'HD2', 'HE1', ) )
        self.add_standard('ILE', ( 'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', 'H', 'HA', 'HB', '1HG1', '2HG1', '1HG2', '2HG2', '3HG2', '1HD1', '2HD1', '3HD1',) )
        self.add_standard('LEU', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'H', 'HA', '1HB ', '2HB ', 'HG', '1HD1', '2HD1', '3HD1', '1HD2', '2HD2', '3HD2',) )
        self.add_standard('LYS', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ', 'H', 'HA', '1HB ', '2HB ', '1HG ', '2HG ', '1HD ', '2HD ', '1HE ', '2HE ', '1HZ ', '2HZ ', '3HZ ', ) )
        self.add_standard('MET', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE', 'H', 'HA', '1HB ', '2HB ', '1HG ', '2HG ', '1HE ', '2HE ', '3HE ',) )
        self.add_standard('PHE', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'H', 'HA', '1HB ', '2HB ', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ', ) )
        self.add_standard('PRO', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', '*NV', 'HA', '1HB ', '2HB ', '1HG ', '2HG ', '1HD ', '2HD ', ) )
        self.add_standard('SER', ( 'N', 'CA', 'C', 'O', 'CB', 'OG', 'H', 'HA', '1HB ', '2HB ', 'HG',) )
        self.add_standard('THR', ( 'N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', 'H', 'HA', 'HB', 'HG1', '1HG2', '2HG2', '3HG2',) )
        self.add_standard('TRP', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'H', 'HA', '1HB ', '2HB ', 'HD1', 'HE1', 'HE3', 'HZ2', 'HZ3', 'HH2', ) )
        self.add_standard('TYR', ( 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH', 'H', 'HA', '1HB ', '2HB ', 'HD1', 'HD2', 'HE1', 'HE2', 'HH', ) )
        self.add_standard('VAL', ( 'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'H', 'HA', 'HB', '1HG1', '2HG1', '3HG1', '1HG2', '2HG2', '3HG2',) )

        self.standard_patch()

    def add_standard(self, name, atoms):
        '''Add a "standard" amino acid residue'''
        self.add_residue(name, atoms)
        self.standard[name] = True

    def add_residue(self, name, atoms):
        '''Add a residue'''
        if name in self.atom_names:
            raise ValueError('Residue type "'+name+'" already exists!')
        self.atom_names[name] = tuple(atoms)

    def standard_patch(self):
        '''Make common patches for standard residues'''
        for base in self.standard: # Assuming they're not already patched
            base_atoms = self.atom_names[base]
            # Nterm
            if 'H' in base_atoms:
                hpos = base_atoms.index('H')
                nterm = base_atoms[:hpos] + ('1H  ','2H  ','3H  ') + base_atoms[hpos+1:]
                self.add_residue(base + '_p:NtermProteinFull', nterm )
            # Cterm
            if 'O' in base_atoms:
                opos = base_atoms.index('O')
                oterm = base_atoms[:opos+1] + (' OXT',) + base_atoms[opos+1:]
                self.add_residue(base + '_p:CtermProteinFull', oterm )
            # Cutpoint lower
            if 'O' in base_atoms:
                opos = base_atoms.index('O')
                lower = base_atoms[:opos+1] + ('*OVL1','*OVL2') + base_atoms[opos+1:]
                self.add_residue(base + '_p:protein_cutpoint_lower', lower )
            # Cutpoint upper
            if 'O' in base_atoms:
                opos = base_atoms.index('O')
                upper = base_atoms[:opos+1] + ('*OVR1',) + base_atoms[opos+1:]
                self.add_residue(base + '_p:protein_cutpoint_upper', upper )

    def get_atoms(self, name):
        '''Get atom names for a given residue type. Translate single letter codes. Return None if type not found.'''
        if len(name) == 1:
            name = self.ONE2THREE.get(name,name)
        if name not in self.atom_names:
            return None
        else:
            return self.atom_names[name]

    def get_resn(self, name):
        '''Get the three-letter residue name for a given input residue name.'''
        if len(name) == 1:
            name = self.ONE2THREE.get(name,name)
        name3 = name.split('_p')[0]
        if len(name3) < 3:
            name3 = "%3s" % name3
        if len(name3) > 3:
            print('WARNING - base typename for "'+name+'" ("'+name3+'") exceeds 3 charachters - truncating.')
            name3 = name3[:3]
        return name3

class PDBFile:
    """A class to hold a PDB file representation."""
    def __init__(self):
        self.lines = []
        self.lastatom = 0
        self.chainendings = []

    def add_atom(self, aname, resn, resi, x, y, z ):
        self.lastatom += 1
        chain = 0
        if self.chainendings:
            for e in self.chainendings:
                if resi > e:
                    chain += 1
                else:
                    break
        chain = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'[chain%52]
        self.lines.append( "ATOM  %5d %4.4s %3.3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00              " % ( self.lastatom, "%-3s"%aname, resn, chain, resi, x, y, z ) )

    def chain_endings(self, endings):
        self.chainendings = [int(e) for e in endings]
        self.chainendings.sort()

    def get_file(self):
        '''Returns the represented PDB file.'''
        return '\n'.join(self.lines)

DATABASE = ResidueDatabase()

class CoordinateDecoder:
    def __init__(self, endianout = None):
        self.endian = endianout

    def decode(self, encoded):
        '''Main entry point for residue coordinate decoding - handles endian detection.'''
        #The remembering of the encoding of the first endianness test is intentional
        if self.endian is not None:
           return self.decode_internal(encoded)
        else:
            self.endian = "<" # Intel little-endian is a good first guess
            coords = self.decode_internal(encoded)
            # Test for big endian encoding.
            if 0.25 <= self.dist2(coords[0],coords[1]) <= 4.0 and 0.25 <= self.dist2(coords[1],coords[2]) <= 4.0:
                return coords
            else:
                #Retry with big-endian
                self.endian = ">"
                coords2 = self.decode_internal(encoded)
                if 0.25 <= self.dist2(coords2[0],coords2[1]) <= 4.0 and 0.25 <= self.dist2(coords2[1],coords2[2]) <= 4.0:
                    print "NOTICE - autodetected a big-endian silent file - autoconverting."
                    return coords2
                else:
                    print "NOTICE - bond lengths are a bit off, but assuming standard format."
                    self.endian = "<" # reset endian
                    return coords

    def decode_internal(self, encoded):
        '''Internal residue coordiante decoding - assumes endian is set.'''
        coords = []
        decoded = [self.decode_char(c) for c in encoded]
        for i in range( int(len(decoded)/16) ):
            coords.append( self.decode_coord( decoded[i*16:(i+1)*16] ) )
        return coords

    def decode_coord(self, decoded):
        '''Decodes a single atom xyz coordinate.'''
        assert len(decoded) == 16
        shifted = self.bitshift(decoded[0:4]) + self.bitshift(decoded[4:8]) + self.bitshift(decoded[8:12]) + self.bitshift(decoded[12:16]) # tuple addition
        assert len(shifted) == 12
        buff = struct.pack(12*'B', *shifted )
        vals = struct.unpack(self.endian + 3 * 'f', buff ) # Always uses IEEE 754 binary32 representation
        assert len(vals) == 3
        return vals

    def decode_char(self, char):
        '''Decodes a single input stream characher'''
        assert len(char) == 1
        if 'A' <= char <= 'Z':
            return ord(char) - ord('A')
        elif 'a' <= char <= 'z':
            return ord(char) - ord('a') + 26
        elif '0' <= char <= '9':
            return ord(char) - ord('0') + 52
        elif '+' == char:
            return 62
        elif '/' == char:
            return 63
        else:
            msg = "Attempted to decode unknown character '"+char+"'"
            print(msg)
            raise ValueError(msg)

    def bitshift(self, charset ):
        assert len(charset) == 4
        total = charset[0] | (charset[1] << 6) | (charset[2] << (2*6)) | (charset[3] << (3*6))
        return ( (total & 0xFF), ((total >> 8) & 0xFF), ((total >> 16) & 0xFF), )

    def dist2(self, a, b):
        assert len(a) == 3 and len(b) == 3
        return (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2

#######

def parse_res(tag, resi, residue, coords, pdb):
    atomnames = DATABASE.get_atoms(residue)
    resn = DATABASE.get_resn(residue)
    if atomnames is None:
        atomnames = DATABASE.get_atoms(resn)
        if atomnames is None:
            print("WARNING - Can't find data for "+residue+" (residue number "+str(resi)+" in "+tag+") - typing as unknown" )
            atomnames = tuple( 'U'+str(i) for i in range(1, len(coords) + 1) )
        else:
            print("WARNING - Can't find data for "+residue+" (residue number "+str(resi)+" in "+tag+") - typing as "+ resn )

    if len(atomnames) > len(coords):
        print("WARNING - Not enough encoded atoms for "+residue+" (residue number "+str(resi)+" in "+tag+") - doing my best (there probably will be mistypes)" )
    if len(atomnames) < len(coords):
        print("WARNING - More encoded atoms for "+residue+" (residue number "+str(resi)+" in "+tag+") than expected - typing extras as unknown (there probably will be mistypes)" )
        atomnames = atomnames + tuple( 'U'+str(i) for i in range(1, len(coords) - len(atomnames) + 1) )

    for i, (x,y,z) in enumerate(coords):
        assert i < len(atomnames)
        name = atomnames[i]
        if name.startswith('*'):
            continue
        pdb.add_atom( name, resn, resi, x, y, z )

def parse_sequence( annotated ):
    parsed = []
    i = 0
    while i < len(annotated):
        name = annotated[i]
	#Sometimes the last amino acid doesn't have additional type specifications.
        if i+1 < len(annotated) and annotated[i+1] == "[":
            i=i+2
            full = []
            while i < len(annotated) and annotated[i] != "]":
                full.append(annotated[i])
                i = i+1
            name = ''.join(full)
        parsed.append(name)
        i = i+1
    return parsed

def load_struct(tag, seqline, f): # f should be sitting at the start of the line after the ANNOTATED_SEQUENCE line.
    pdb = PDBFile()
    decoder = CoordinateDecoder() # A new decoder for each structure - as a file might be a concatenation of multiple files with different endianness
    resis = parse_sequence(seqline.split()[1])
    curres = 1
    while True:
        line = f.readline()
        if not line:
            break # End of file
        line = line.strip()
        if not line:
            continue # blank line

        if not line or (line.startswith("REMARK") or
                        line.startswith("SCORE: ") or
                        line.startswith("FOLD_TREE ") or
                        line.startswith("RT") or
                        line.startswith("SEQUENCE:") or
                        line.startswith("JUMP") or
                        line.startswith("SYMMETRY_INFO") ):
             continue
        if line.startswith("CHAIN_ENDINGS"):
            line = line.split()
            if line[-1] != tag:
                print("WARNING - tag mismatch in coordinates. Possible missing header " + tag + " and " + line[-1])
                return
            pdb.chain_endings(line[1:-1])
            continue
        if line.startswith("ANNOTATED_SEQUENCE:"):
            print("WARNING - duplicated structure lines - skipping " + tag )
            return
        #If we get here, it should be a coordinate line.
        line = line.split()
        if len(line) != 2:
            print("WARNING - bad format in coordinate line - skipping " + tag )
            return
        if line[-1] != tag:
            print("WARNING - tag missmatch in coordinates. Possible missing header " + tag + " and " + line[-1])
            return
        ss, coords = line[0][0], decoder.decode(line[0][1:])
        parse_res(tag, curres, resis[curres-1], coords, pdb)
        if curres == len(resis):
            break
        curres += 1

    return pdb.get_file()

def checktags(tag, rextags):
    '''Return true if tag matches a (compiled regex) entry in rextags, false if not.'''
    if len(rextags) == 0 :
        return True
    for rex in rextags:
        if rex.match(tag):
            return True
    return False

def silent_load(filename, *tags, **kwds): # tags = empty loads all the strucutres. (kwds needed to eat a superflous _self parameter)
    """silent_load filename, [tags]

    Load the given tags from the given Rosetta binary-format silent file into appropriately named structures in Pymol.
    Tags can be zero or more comma separated tag. If no tags are given, load all availible tags.
    Python-style regex expansion can be used to match multiple tags. Regexes must match the entire tag to be accepted.
    """
    # Add start and end so we don't get partial matches.
    # Strip out whitespace and quote marks (shouldn't occur in tags or regex, but might occur due to terminal quoting)
    rextags = [ re.compile( "^" + t.strip().strip('"'+"'") + "$") for t in tags ]
    with open(filename) as f:
        while True:
            line = f.readline()
            if not line:
                break  # end of file
            line = line.strip()
            if not line:
                continue
            #In practice, ANNOTATED_SEQEUENCE always comes before coordiantes
            if line.startswith("ANNOTATED_SEQUENCE:"):
                tag = line.split()[-1]
                if checktags(tag, rextags):
                    PDBLINES = load_struct(tag, line, f)
                    if PDBLINES:
                        if cmd is not None:
                            # Load the PDB into Pymol
                            print("-- Loading structure " + tag)
                            cmd.read_pdbstr( PDBLINES, tag )
                        else:
                            outfilename = tag+".pdb"
                            print("-- Saving file " + outfilename)
                            with open(outfilename,'w') as of:
                                of.write(PDBLINES)
                            pass

def silent_tags(filename):
    """silent_tags filename

    Show all the availible tags associated with the given Rosetta binary-format silent file."""
    tags = []
    with open(filename) as f:
        for line in f:
            if line.startswith("ANNOTATED_SEQUENCE:"):
                tag = line.split()[-1]
                tags.append(tag)
    tags.sort()
    width = max([len(t) for t in tags]) + 4
    number = int(80/width)
    for i, t in enumerate(tags):
        sys.stdout.write( ("%-"+str(width)+"s") % t )
        if i % number == number - 1:
            sys.stdout.write('\n')
    sys.stdout.write('\n')


if cmd is not None:
    # Running under pymol
    cmd.extend("silent_load",silent_load)
    cmd.extend("silent_tags",silent_tags)
elif __name__ == "__main__":
    # Running as a script
    silent_load(sys.argv[1], *sys.argv[2:])
else:
    pass
