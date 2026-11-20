#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

'''Utilities for reading and manipulating MDL Molfiles (.mol, .mdl) and SD files (.sdf)
as well as Tripos MOL2 (.mol2) files.

Molfiles are common, simple, text-based listings of atom positions, connectivity, and bond order.
MDL molfiles specify a single conformer of a single molecule.
SD files are basically multiple molfiles separated by "$$$$".
See the official specification at http://www.mdl.com/downloads/public/ctfile/ctfile.jsp.

MOL2 files are very similar in spirit to .mol/.sdf, but are not column-oriented (which is nice).
MOL2 files originated with Tripos (a company) and Sybyl (software).
See the official specification at http://tripos.com/data/support/mol2.pdf

Author: Ian W. Davis
'''

from __future__ import print_function
import os, sys, math, copy, gzip

try: set
except: from sets import Set as set

#open file normally or with gzip depending on the file extension
def gz_open(file,mode):
    filename, extension = os.path.splitext( file )
    if extension == ".gz":
        return gzip.open(file,mode)
    else:
        return open(file,mode)

'''Dictionary of Rosetta atom names to plausible PDB atom names (element only)'''
rosetta_to_pdb_names = {'CNH2':' C  ', 'COO ':' C  ', 'CH1 ':' C  ',
    'CH2 ':' C  ', 'CH3 ':' C  ', 'aroC':' C  ', 'Ntrp':' N  ', 'Nhis':' N  ',
    'NH2O':' N  ', 'Nlys':' N  ', 'Narg':' N  ', 'Npro':' N  ', 'OH  ':' O  ',
    'ONH2':' O  ', 'OOC ':' O  ', 'S   ':' S  ', 'Nbb ':' N  ', 'CAbb':' C  ',
    'CObb':' C  ', 'OCbb':' O  ', 'Hpol':' H  ', 'Hapo':' H  ', 'Haro':' H  ',
    'HNbb':' H  ', 'Phos':' P  '}

class Molfile:
    # title         molecule name or identifier
    # comment       optional comment describing the molecule
    # atoms[]       list of Atom objects, in file order
    # bonds[]       list of Bond objects, in file order
    # footer[]      zero or more footer lines, like M____ records (with newlines)

    def __init__(self):
        self.title = ""
        self.comment = ""
        self.atoms  = []
        self.bonds  = []
        self.footer = []

class Atom:
    # x, y, z           coordinates
    # name              atom name without whitespace
    # sybyl_type        one of the SYBYL atom type strings, default "Any"
    # rsd_id            residue or substructure ID number, default 1
    # rsd_name          residue or substructure name, default "RSD"
    # elem              element, 1-2 letters upper case
    # is_H              is this atom a hydrogen?
    # bonds[]           bonds for this atom, with self == a1
    # heavy_bonds[]     ditto, but only bonds TO non-H atoms
    # is_ring           is this atom part of a ring/cycle?
    # ring_size         smallest ring this atom is in, or zero
    # partial_charge    may be None, else a real number

    def __init__(self, x, y, z, name, elem=None):
        self.x = x
        self.y = y
        self.z = z
        self.name = name
        self.sybyl_type = "Any"
        self.rsd_id = 1
        self.rsd_name = "RSD"
        key = "%-4s" % name # left-justified like PDB name
        if elem is not None:
            self.elem = elem
        elif key in rosetta_to_pdb_names:
            self.elem = rosetta_to_pdb_names[key].strip()
        elif name[0].isalpha():
            if len(name) >= 2 and name[1].isalpha():
                self.elem = name[0:2]
            else:
                self.elem = name[0]
        else:
            self.elem = name[1]
        self.elem = self.elem.upper()
        self.is_H = (self.elem == "H")
        self.bonds = []
        self.heavy_bonds = []
        self.is_ring = False
        self.ring_size = 0
        self.partial_charge = None
        self.formal_charge = 0

    def copy(self):
        '''Return a semi-shallow copy of this Atom, with bonds[] and heavy_bonds[] set to empty lists.'''
        dup = copy.copy(self)
        for name, val in dup.__dict__.iteritems():
            dup.__dict__[name] = copy.copy(val)
        dup.bonds = []
        dup.heavy_bonds = []
        return dup

    def __str__(self):
        return "(%4s %2s %8.3f %8.3f %8.3f)" % (self.name, self.elem, self.x, self.y, self.z)

    def __lt__(self, other):
        '''A quick-and-dirty sorting method, to allow for decorate-sort-undecorate sorting of Atom lists.'''
        return self.elem < other.elem

def index_atoms(atoms):
    '''Returns a dictionary of the 1-based indices of given Atoms.'''
    ai = {}
    for i,a in enumerate(atoms):
        ai[a] = i+1
    return ai

class Bond:
    # static constants:
    SINGLE      = 1
    DOUBLE      = 2
    TRIPLE      = 3
    AROMATIC    = 4 # loosely speaking -- potentially includes amide and conjugated bonds

    # a1, a2        the bonded atoms
    # order         single, double, triple, aromatic
    # is_ring       is this bond part of a ring/cycle?
    # ring_size     smallest ring this atom is in, or zero
    # mirror        same bond with a1 and a2 swapped

    def __init__(self, a1, a2, order, mirror=None):
        self.a1 = a1
        self.a2 = a2
        self.order = order
        self.is_ring = False
        self.ring_size = 0
        if mirror is None:
            self.mirror = Bond(a2, a1, order, self)
            # Go ahead and add to atoms' bond lists so we don't forget to do it
            a1.bonds.append(self)
            if not a2.is_H: a1.heavy_bonds.append(self)
            # Use mirror for a2 so it will be in expected "a1" position
            a2.bonds.append(self.mirror)
            if not a1.is_H: a2.heavy_bonds.append(self.mirror)
        elif self.a1 == mirror.a2 and self.a2 == mirror.a1 and self.order ==  mirror.order:
            self.mirror = mirror
        else:
            raise ValueError("bad mirror")

    def __str__(self):
        if   self.order == Bond.SINGLE:   bnd = " --- "
        elif self.order == Bond.DOUBLE:   bnd = " === "
        elif self.order == Bond.TRIPLE:   bnd = " ### "
        elif self.order == Bond.AROMATIC: bnd = " ::: "
        else: bnd = " ??? "
        s = str(self.a1) + bnd + str(self.a2)
        if self.is_ring: s += " RING %i" % self.ring_size
        return s

def find_rings(bonds):
    '''A bond is in a ring iff there is an alternate path from a1 to a2.

    An atom is in a ring iff it participates in any ring bonds.
    Calling this function sets is_ring for all atoms and bonds.'''
    for bond in bonds:
        # For each bond, do a breadth-first search of the bond graph,
        # not using the current bond under consideration.
        visited = set()
        visited.add(bond.a1)
        to_visit = []
        min_dist = {}
        # Start with all neighbors except a2
        for bnd in bond.a1.bonds:
            if bnd.a2 != bond.a2:
                to_visit.append(bnd.a2)
                min_dist[ bnd.a2 ] = 1
        # Do breadth-first search
        while len(to_visit) > 0:
            a = to_visit.pop(0)
            visited.add(a)
            for bnd in a.bonds:
                if bnd.a2 not in visited:
                    to_visit.append(bnd.a2)
                    min_dist[ bnd.a2 ] = min_dist[a] + 1 # guaranteed to be min b/c search is breadth-first
        # Did we reach the other side?
        def posmin(a, b):
            if a > 0:
                if b > 0: return min(a, b)
                else: return a
            elif b > 0: return b
            else: assert False, "Expected at least one positive argument"
        if bond.a2 in visited:
            ring_size = min_dist[ bond.a2 ] + 1 # one more bond to close the ring
            bond.is_ring = True
            bond.ring_size = posmin(ring_size, bond.ring_size)
            bond.mirror.is_ring = True
            bond.mirror.ring_size = bond.ring_size
            bond.a1.is_ring = True
            bond.a1.ring_size = posmin(ring_size, bond.a1.ring_size)
            bond.a2.is_ring = True
            bond.a2.ring_size = posmin(ring_size, bond.a2.ring_size)

def file_or_filename(func):
    '''A decorator for functions that interchangably take a file or filename as their first argument.'''
    def g(f, *args, **kwargs):
        if isinstance(f, str):
            f = gz_open(f, 'r')
            ret = func(f, *args, **kwargs)
            #f.close()
            # Can't close files when we do it like this,
            # because the generator may not be used until much later.
            # Just have to wait for it to go out of scope and get GC'd.
            return ret
        else:
            return func(f, *args, **kwargs)
    # Preserve some signature info from wrapped function
    g.__name__ = func.__name__
    g.__doc__ = func.__doc__
    g.__dict__.update(func.__dict__)
    return g

def assign_mdl_charges(line, atoms):
    '''Assigns formal charges based on the MDL M  CHG line.'''
    nentries = int(line[6:9])
    line = line[9:].split()
    if( nentries*2 != len(line) ):
        print( "Warning: Malformed charge line. Ignoring." )
        return
    for i in range(nentries):
        atom = int(line[2*i]);
        charge = int(line[2*i + 1])
        atoms[atom-1].formal_charge = charge # 1 based to zero-based indexing

def read_mdl_molfile(f, do_find_rings=True):
    '''Reads a molfile and returns a Molfile object.

    f may be a file name or file handle.'''
    if isinstance(f, str):
        f = gz_open(f, 'rU')
        ret = read_mdl_molfile(f)
        f.close()
        return ret
    molfile = Molfile()
    atoms = molfile.atoms
    bonds = molfile.bonds
    # Discard first three lines, which are comments
    molfile.title = (f.readline().strip() + " " + f.readline().strip()).strip()
    molfile.comment = f.readline().strip()
    # Next line defines number of atoms, number of bonds, and format version
    fmt = f.readline()
    if fmt == "": return None # EOF check for reading reading SDF entries
    # Concession to sloppy files that are missing the V2000 signature ... sigh.
    if len(fmt) >= 39 and fmt[34:39] != "V2000":
        raise ValueError("can only read V2000 format files")
    n_atoms = int(fmt[0:3])
    n_bonds = int(fmt[3:6])
    # Read prescribed number of atom lines
    for i in range(n_atoms):
        line = f.readline()
        if line == "": raise ValueError("premature end of file: atom %i" % i)
        atoms.append( Atom(
            float(line[ 0:10]),
            float(line[10:20]),
            float(line[20:30]),
            # technically it's just 3 chars, but we often need/use 4 chars
            line[30:35].strip() ))
    # Read prescribed number of bond lines
    for i in range(n_bonds):
        line = f.readline()
        if line == "": raise ValueError("premature end of file: bond %i" % i)
        a1 = int(line[0:3]) - 1
        a2 = int(line[3:6]) - 1
        order = int(line[6:9])
        bond = Bond(atoms[a1], atoms[a2], order)
        bonds.append(bond)
    # Read any remaining lines to save for later
    while True:
        line = f.readline()
        if line == "" or line.startswith("$$$$"): break
        elif line.startswith("M  END"): continue
        elif line.startswith("M  CHG"):
            assign_mdl_charges(line, atoms)
        else:
            molfile.footer.append(line)
    if do_find_rings: find_rings(bonds)
    return molfile

def write_mdl_charges(f, molfile, atom_idx):
    charged_atoms = []
    for a in molfile.atoms:
        if a.formal_charge != 0:
            charged_atoms.append(a)
    # CHG lines can only have a maximum of 8 atoms
    b = 0
    e = 8
    while b < len(charged_atoms):
        end = min(e, len(charged_atoms) )
        nentries = end-b
        assert nentries <= 8
        f.write("M  CHG%3d" % nentries)
        for n in range(b,end):
            f.write(" %3d %3d" % (atom_idx[charged_atoms[n]], charged_atoms[n].formal_charge) )
        f.write("\n")
        b += 8
        e += 8

def write_mdl_molfile(f, molfile):
    '''Writes a Molfile object to a file.

    If atoms and/or bonds have been re-ordered since read-in,
    this code will renumber everything appropriately.

    f may be a file name or file handle.'''
    if isinstance(f, str):
        f = gz_open(f, 'w')
        write_mdl_molfile(f, molfile)
        f.close()
        return
    f.write(molfile.title+"\n")
    f.write("\n")
    f.write(molfile.comment+"\n")
    # Overwrite atom and bond counts, just in case they changed
    fmt = "%3i%3i  0     1  0  0  0  0  0999 V2000\n" % (len(molfile.atoms), len(molfile.bonds))
    f.write(fmt)
    atom_idx = index_atoms(molfile.atoms) # so bonds can look it up later!
    for a in molfile.atoms:
        # For now we discard the rest of the fields, so just make them zeros.
        # Atom name should be three char, but we will often have four.
        # This code matches the "Rosetta standard" but the zeros are misplaced
        # by one column now; may require lenience from other readers
        # but won't introduce spurious differences in our existing molfiles.
        f.write("%10.4f%10.4f%10.4f %-4s 0  0  0  0  0  0  0  0  0  0  0  0\n" % (a.x, a.y, a.z, a.name))
    for b in molfile.bonds:
        # For now we discard the rest of the fields, so just make them zeros
        f.write("%3i%3i%3i  0  0  0  0\n" % (atom_idx[b.a1], atom_idx[b.a2], b.order))
    write_mdl_charges(f, molfile, atom_idx)
    f.writelines(molfile.footer)
    f.write("M  END\n")
    f.flush() # close() ruins StringIO objects!

@file_or_filename
def read_mdl_sdf(f, do_find_rings=True):
    '''Reads an sdf and returns a list of Molfile objects.

    f may be a file name or file handle.'''
    if isinstance(f, str):
        f = gz_open(f, 'rU')
        ret = read_mdl_sdf(f, do_find_rings)
        f.close()
        return ret
    molfiles = []
    while True:
        molfile = read_mdl_molfile(f, do_find_rings)
        if molfile is None: break
        molfiles.append(molfile)
    if len(molfiles) == 0: raise ValueError("no entries in SDF file")
    return molfiles

def write_mdl_sdf(f, molfiles):
    '''Writes a list of Molfile objects to a file.

    If atoms and/or bonds have been re-ordered since read-in,
    this code will renumber everything appropriately.

    f may be a file name or file handle.'''
    if isinstance(f, str):
        f = gz_open(f, 'w')
        write_mdl_sdf(f, molfiles)
        f.close()
        return
    for molfile in molfiles:
        write_mdl_molfile(f, molfile)
        f.write("\n$$$$\n")
    f.flush()

@file_or_filename
def read_tripos_mol2(f, do_find_rings=True):
    '''Reads a mol2 and returns a list of Molfile objects.

    f may be a file name or file handle.'''
    if isinstance(f, str):
        f = gz_open(f, 'rU')
        ret = read_tripos_mol2(f, do_find_rings)
        f.close()
        return ret
    # Line generator function to deal with comments and continuations
    line_num = [0] # stupid Python can't assign to outer variables!
    def read_mol2_lines(f):
        full_line = ""
        while True:
            line = f.readline()
            line_num[0] += 1 # stupid Python can't assign to outer variables!
            if line == "":
                assert full_line == "" # can't continue into EOF
                return # EOF
            elif line.startswith("#"):
                assert full_line == "" # can't have comments inside continuations
                continue # skip comments
            else:
                if full_line == "": full_line = line.rstrip()
                else: full_line = (full_line + " " + line).rstrip()
                if full_line.endswith("\\"): full_line = full_line[:-1] # remove trailing backslash
                elif full_line == "": pass # skip empty lines
                else:
                    yield full_line
                    full_line = ""
    line_itr = read_mol2_lines(f)
    molfiles = []
    mode = ""
    for line in line_itr:
        if line == "@<TRIPOS>MOLECULE":
            molfile = Molfile()
            molfiles.append(molfile)
            mode = line
            linecnt = 0
        elif line == "@<TRIPOS>ATOM":
            atom_indices = {}
            mode = line
        elif line == "@<TRIPOS>BOND":
            mode = line
        elif line.startswith("@"):
            mode = line
            molfile.footer.append(line+"\n")
            molfile.footer.append(line[1:]+"\n")
        elif mode == "@<TRIPOS>MOLECULE":
            linecnt += 1
            if linecnt == 1: molfile.title = line.strip()
            elif linecnt == 6: molfile.comment = line.strip()
        elif mode == "@<TRIPOS>ATOM":
            f = line.split()
            assert len(f) >= 6, "Missing fields on line %i" % line_num[0]
            idx = int(f[0])
            name = f[1]
            x = float(f[2])
            y = float(f[3])
            z = float(f[4])
            elem = f[5].split(".")[0]
            atom = Atom(x, y, z, name, elem)
            atom.sybyl_type = f[5]
            if len(f) >= 7:
                try:
                    atom.rsd_id = int(f[6])
                except ValueError:
                    print("Pymol produced mol2 files don't conform to the standard - ignore its error gracefully")
                    f = f[:6]
            if len(f) >= 8:
                atom.rsd_name = f[7]
            if len(f) >= 9:
                atom.partial_charge = float(f[8])
            atom_indices[idx] = atom
            molfile.atoms.append(atom)
        elif mode == "@<TRIPOS>BOND":
            f = line.split()
            assert len(f) >= 4, "Missing fields on line %i" % line_num[0]
            atom1 = atom_indices[int(f[1])]
            atom2 = atom_indices[int(f[2])]
            if f[3] == "1": order = Bond.SINGLE
            elif f[3] == "2": order = Bond.DOUBLE
            elif f[3] == "3": order = Bond.TRIPLE
            elif f[3] == "ar" or f[3] == "am": order = Bond.AROMATIC
            elif f[3] == "du" or f[3] == "un" or f[3] == "nc":
                print( "NOTE: Bond of order '%s' treated as single bond." % f[3] )
                order = Bond.SINGLE
            else:
                raise ValueError("Unrecognized bond order '%s' on line %i" % ( f[3], line_num[0] ))
            bond = Bond(atom1, atom2, order)
            molfile.bonds.append(bond)
        else:
            molfile.footer.append(line+"\n")
    if len(molfiles) == 0: raise ValueError("no entries in MOL2 file")
    if do_find_rings:
        for molfile in molfiles: find_rings(molfile.bonds)
    return molfiles

def write_tripos_mol2(f, molfiles):
    '''Writes a list of Molfile objects to a file.

    If atoms and/or bonds have been re-ordered since read-in,
    this code will renumber everything appropriately.

    f may be a file name or file handle.

    This function doesn't preserve everything, notably substructure records and amide bond types.'''
    if isinstance(f, str):
        f = gz_open(f, 'w')
        write_tripos_mol2(f, molfiles)
        f.close()
        return
    for molfile in molfiles:
        f.write("@<TRIPOS>MOLECULE\n")
        f.write(molfile.title+"\n")
        unique_rsd_ids = set([ (a.rsd_id, a.rsd_name) for a in molfile.atoms ])
        f.write("%6i %6i %6i\n" % (len(molfile.atoms), len(molfile.bonds), len(unique_rsd_ids)))
        f.write("SMALL\n")
        if len([a for a in molfile.atoms if a.partial_charge is not None]) > 0:
            f.write("USER_CHARGES\n")
        else:
            f.write("NO_CHARGES\n")
        f.write("\n") # status flags
        f.write(molfile.comment+"\n")
        f.write("@<TRIPOS>ATOM\n")
        atom_idx = index_atoms(molfile.atoms) # so bonds can look it up later!
        for a in molfile.atoms:
            f.write("%8i %-8s %9.4f %9.4f %9.4f %-5s %5i %-8s" % (atom_idx[a], a.name, a.x, a.y, a.z, a.sybyl_type, a.rsd_id, a.rsd_name))
            if a.partial_charge is not None: f.write(" %9.4f" % a.partial_charge)
            f.write("\n")
        f.write("@<TRIPOS>BOND\n")
        for i, b in enumerate(molfile.bonds):
            f.write("%8i %8i %8i " % (i+1, atom_idx[b.a1], atom_idx[b.a2]))
            if b.order == Bond.AROMATIC: f.write("ar")
                # Can you ever have an amide N with two aromatic bonds to it?  Maybe in a ring?  Then this would fail...
                #if b.a1.sybyl_type == "N.am" or b.a2.sybyl_type == "N.am": f.write("am")
                #else: f.write("ar")
            else: f.write("%i" % b.order)
            f.write("\n")
        f.writelines(molfile.footer)
    f.flush() # close() ruins StringIO objects!

def sort_for_rosetta(molfile):
    '''Sorts the atoms and bonds of a Molfile into the order preferred by Rosetta:
    Heavy atoms precede hydrogens, and bonds define a tree.'''
    m = molfile
    # Sort atoms so all H go to the end
    m.atoms.sort(lambda a,b: cmp(a.is_H, b.is_H))
    # Direct bonds so lesser atom index comes first
    ai = index_atoms(m.atoms) # atom indices
    bs = [];
    for b in m.bonds:
        if(ai[b.a1] <= ai[b.a2]): bs.append(b)
        else: bs.append(b.mirror)
    # Sort bonds by first atom index, then by second atom index
    # This is a quality Rosetta (used to) expect:
    # (1) the first two atoms are heavy, and
    # (2) every bond goes from a known atom to a (new or known) atom.
    # This defines an implicit tree.
    def bond_cmp(b1, b2):
        c = cmp(ai[b1.a1], ai[b2.a1])
        if c == 0: c = cmp(ai[b1.a2], ai[b2.a2])
        return c
    bs.sort(bond_cmp)
    m.bonds = bs

def strip_H(molfile, pred=lambda x: x.is_H):
    '''Removes hydrogen atoms from a Molfile object, modifying it in place.
    Can remove an arbitrary set of atoms by also passing a predicate
    that takes an Atom and returns True for atoms to be removed.'''
    def strip_bonds(bonds):
        return [b for b in bonds if not pred(b.a1) and not pred(b.a2)]
    molfile.atoms = [a for a in molfile.atoms if not pred(a)]
    for a in molfile.atoms:
        a.bonds = strip_bonds(a.bonds)
        a.heavy_bonds = strip_bonds(a.heavy_bonds)
    molfile.bonds = strip_bonds(molfile.bonds)

def pdb_pad_atom_name(atom):
    '''Returns the atom name padded with whitespace to match the PDB conventions.'''
    if len(atom.elem) == 1 and len(atom.name) <= 3: return " %-3s" % atom.name
    else: return "%-4s" % atom.name

def uniquify_atom_names(atoms, force=False):
    '''If force is true, rename/number all atoms.
    Otherwise, only rename/number atoms with non-unique names.
    Return True iff atoms are renamed.'''
    all_atom_names = [a.name.strip() for a in atoms]
    atom_names = set(all_atom_names)
    if force:
        dup_names = atom_names
        atom_names = set()
    elif len(atom_names) == len(atoms):
        return False # all names are unique, no changes
    else:
        unique_names = set()
        dup_names = set()
        for name in all_atom_names:
            if name not in unique_names: unique_names.add(name)
            else: dup_names.add(name)
        del unique_names # not accurately named - each atom name in here once
        atom_names -= dup_names
    for atom in atoms:
        if atom.name.strip() not in dup_names: continue
        i = 1
        while True:
            name = "%s%i" % (atom.elem[0:2], i)
            if name not in atom_names: break
            i += 1
        atom.name = name
        atom_names.add(name)
    return True # names modified

#def main(argv):
    #ms = list(read_tripos_mol2("1aq1.mol2"))
    #write_mdl_sdf("tmp.sdf", ms)
    #m = read_mdl_molfile("1aq1_ligand.mol")
    #for b in m.bonds: print b
    #make_chis_and_icoor(m.atoms, m.bonds)
    #m = read_mdl_molfile("1aq1_ligand.mol")
    #sort_for_rosetta(m)
    #write_mdl_molfile("tmp.mol", m)

#if __name__ == "__main__":
    #sys.exit(main(sys.argv[1:]))
