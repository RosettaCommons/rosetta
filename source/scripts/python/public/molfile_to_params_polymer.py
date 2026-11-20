#!/usr/bin/env python
'''
Functions and executable for taking a small molecule from an MDL Molfile
and splitting it into one or more .params files for Minirosetta.
See main() for usage or run with --help.

Author: Ian W. Davis, Oanh Vu, Benjamin P. Brown
'''

from __future__ import print_function

import os, sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

# Magic spell to make sure Rosetta python libs are on the PYTHONPATH:
sys.path.append( os.path.abspath( sys.path[0] ) )

from rosetta_py.io.mdl_molfile import *
#from rosetta_py.utility.rankorder import argmin
from rosetta_py.utility import r3
import argparse

# import helper functions
from molfile_to_params_polymer.atom_functions import *
from molfile_to_params_polymer.bond_functions import *
from molfile_to_params_polymer.fragment_functions import *
from molfile_to_params_polymer.polymer_functions import *
from molfile_to_params_polymer.IO_functions import *

# Features from Python 2.5 that we want to use:
if not hasattr(__builtins__, "any"):
    def any(itr):
        for el in itr:
            if el: return True
        return False

if not hasattr(__builtins__, "all"):
    def all(itr):
        for el in itr:
            if not el: return False
        return True

def main(argv):
    """
Converts a small molecule in an MDL Molfile with "M SPLT" and "M ROOT"
records into a series of .params residue definition files for Rosetta.
Also writes out the ligand conformation as PDB HETATMs.
If an SD file is given as input instead, the first entry is used for
generating topology / parameter files, and they all are used for
generating PDB-style coordinates in separate, numbered files.
Multiple models may also be supplied in MOL2 format, which does not support
M ROOT and M SPLT records but does allow for partial charges.
File type is deduced from the extension.

To divide a ligand into fragments by "breaking" bonds (optional):
M SPLT atom_no1 atom_no2

To specify a root atom for a ligand fragment (optional):
M ROOT atom_no

Expects that the input ligand has already had aromaticity "perceived",
i.e. that it contains aromatic bonds rather than alternating single and double
bonds (Kekule structure).

Optionally writes a kinemage graphics visualization of the atom tree,
neighbor atom selection, fragments, etc -- very helpful for debugging
and for visualizing exactly what was done to the ligand.
    """

    parser = argparse.ArgumentParser(description='This script creates Rosetta params files'
                                     'from the molecule\'s mol/sdf file')
    parser.add_argument("-i","--input", required=True,
                        help="name of the intput mol/sdf file")

    parser.add_argument("-n", "--name", default="LG",
                        help="name ligand residues NM1,NM2,... instead of LG1,LG2,...",
                        metavar="NM"
                       )
    parser.add_argument("-p", "--pdb",
                        default="LG", # same as --name, see below
                        help="prefix for PDB file names",
                        metavar="FILE"
                       )
    parser.add_argument("-c", "--centroid",
                        default = None, # same as --name, see below,
                        nargs = '+',
                        help = "translate output PDB coords to have given heavy-atom centroid coods",
                        metavar="X,Y,Z"
                       )
    parser.add_argument("-m", "--max-confs", type = int,
                        default=5000, # 400 (default Omega max) * 9 (one sp3 H with -ex1) = 3600
                        help="don't expand proton chis if above this many total confs",
                        metavar="MAX"
                       )
    parser.add_argument("-k", "--kinemage",
                        default = None,
                        help="write ligand topology to FILE",
                        metavar="FILE"
                       )
    parser.add_argument("--clobber",
                        default=False,
                        action="store_true",
                        help="overwrite existing files"
                       )
    parser.add_argument("--no-param",
                        default = False,
                        action = "store_true",
                        help = "skip writing .params files (for debugging)"
                       )
    parser.add_argument("--no-pdb",
                        default = False,
                        action = "store_true",
                        help = "skip writing .pdb files (for if not using PDB rotamers)"
                       )
    parser.add_argument("--all-in-one-pdb",
                        default = False,
                        action = "store_true",
                        help = "writing all pdb files into 1 file (for debugging or PDB rotamers)"
                       )
    parser.add_argument("--use-pdb-rotamers",
                        default = False,
                        action = "store_true",
                        help = "Append PDB_ROTAMERS specification and corresponding file to params file"
                       )
    parser.add_argument("--use-parent-rotamers",
                        default=None,
                        help="One-or three-letter code for the parent amino acid whose backbone and sidechain rotamers "
                       "and rama prepro terms will be used, e.g., 'TRP' or 'W' for tryptophan. All capital letters required. "
                       "No effect for non-polymer params files. Must be one of the 20 canonical amino acids",
                        metavar="PARENT_CAA"
                       )
    parser.add_argument("--polymer",
                        default = True,
                        action = "store_true",
                        help = "write a polymer style param file instead of a ligand param file"
                       )
    parser.add_argument("--peptoid",
                        default = False,
                        action = "store_true",
                        help = "modifier for the polymer flag, " +
                        "adjusts PDB style naming to be correct for peptoids"
                       )
    parser.add_argument("--keep-names",
                        default=False,
                        action="store_true",
                        help="leaves atom names untouched except for duplications"
    )
    parser.add_argument('--partial_charges',
                        default=None,
                        help="file that contains the partial charges of each atom in the input file",
                        metavar="FILE")

    args = parser.parse_args()

    if args.pdb == 'LG': args.pdb = args.name


    # There's a very strong order dependence to these function calls:
    # many depend on atom/bond variables set by earlier calls.
    if args.input.endswith(".mol2"):
        molfiles = read_tripos_mol2(args.input)
    elif args.input.endswith(".mol") or args.input.endswith(".mdl") or args.input.endswith(".sdf"):
        molfiles = read_mdl_sdf(args.input)
    else:
        print( "Unrecognized file type, must be .mol/.sdf or .mol2!" )
        return 6

    m = molfiles[0]
    #m = molfiles
    #print(m)

    ctr = r3.Triple()
    if args.centroid != None:
        if len(args.centroid) == 3:
            ctr = r3.Triple( args.centroid[0], args.centroid[1], args.centroid[2] )
    else:
        print("No valid centroid coords are provided, calculate centroid manually")
        ctr = r3.centroid([a for a in m.atoms if not a.is_H])
        print("Centering ligands at %s" % ctr)
        # If -centroid not given, default is to center like first entry

    add_fields_to_atoms(m.atoms)
    add_fields_to_bonds(m.bonds)
    find_virtual_atoms(m.atoms)
    uniquify_atom_names(m.atoms) # renumber and resname the atoms if there is duplicate atom names
    check_bond_count(m.atoms)
    check_aromaticity(m.bonds)

    if args.polymer:
        polymer_assign_backbone_atom_types(m) # assign backbone atoms
        polymer_assign_ignored_atoms_bonds(m) # ignore the dipeptide atoms
    assign_rosetta_types(m.atoms) # assign rosetta atom names
    assign_mm_types(m.atoms, args.peptoid) # assign mm atom names
    # assign the formal charge based on the footer in the input file
    net_charge = 0.0
    for line in m.footer:
        # KWK's convention
        if line.startswith("M CHG"): net_charge = float(line[5:])
        # Official MDL format (integer only)
        # Oops, this isn't right!  MDL records N atom1 charge1 atom2 charge2 ...
        elif line.startswith("M  CHG"):
            charge_fields = line.split()[3:]
            net_charge = sum(int(c) for i,c in enumerate(charge_fields) if i%2 == 1)
        elif line.startswith("M  POLY_CHG"): net_charge = float(line.split()[-1])
    # If the partial charge info file is provided, then assign those charge values to atoms
    if args.partial_charges != None:
        assign_partial_charges_from_values(m,
                                           read_parital_charge_input(args.partial_charges),
                                           net_charge
                                          )

    # TODO: merge assign_partial_charges_from_values and assign_partial_charges
    assign_partial_charges(m.atoms, args.partial_charges, net_charge) # assign the partial charges to atoms
    assign_rotatable_bonds(m.bonds)
    assign_rigid_ids(m.atoms)
    num_frags = fragment_ligand(m)
    build_fragment_trees(m)
    if args.polymer:
        print("Preforming polymer modifications")
        if args.keep_names:
            for a in m.atoms:
                a.pdb_name = a.name
        else:
            polymer_assign_pdb_like_atom_names_to_sidechain( m.atoms, m.bonds, args.peptoid )
            polymer_assign_backbone_atom_names( m.atoms, m.bonds, args.peptoid )
        # Connections get renamed even if we don't rename other atoms
        polymer_assign_connection_atom_names( m.atoms, m.bonds, args.peptoid )
        # if instructed to write out all pdb files
        if not args.no_pdb :
            for molfile in molfiles[1:]:
                copy_atom_and_bond_info(m, molfile)
                polymer_reorder_atoms(molfile)
        polymer_reorder_atoms(m)
    #uniquify_atom_names(m.atoms)
    assign_internal_coords(m)
    if not args.no_param:
        for i in range(num_frags):
            if num_frags == 1: param_file = "%s.params" % args.name
            else: param_file = "%s%i.params" % (args.pdb, i+1)
            if not args.clobber and os.path.exists(param_file):
                print( "File %s already exists -- aborting!" % param_file )
                print( "Use --clobber to overwrite existing files." )
                return 2
            else:
                if args.polymer:
                    write_poly_param_file(param_file, m, args.name, 1, args.peptoid, args.use_parent_rotamers)
                    print( "Wrote polymer params file %s" % param_file )
                else:
                    write_param_file(param_file, m, args.name, i+1, len(molfiles), args.max_confs)
                    print( "Wrote params file %s" % param_file )
    if args.kinemage is not None:
        if not args.clobber and os.path.exists(args.kinemage):
            print( "File %s already exists -- aborting!" % args.kinemage )
            print( "Use --clobber to overwrite existing files." )
            return 3
        else:
            write_ligand_kinemage(args.kinemage, m)
            print( "Wrote kinemage file %s" % args.kinemage )

    # Info for user
    if args.use_parent_rotamers is not None:
        print("Using parent rotamers: " + str(args.use_parent_rotamers))

    # write out the pdb files
    if not args.no_pdb:
        # check if the same output pdb already exists for writing to 1 pdb
        # if yes, then nothing else to do
        if args.all_in_one_pdb:
            pdb_file_name = "%s_rotamer.pdb" % (args.pdb)
            if not args.clobber and os.path.exists(pdb_file_name):
                print( "File %s already exists -- skip!" % pdb_file_name )
                print( "Use --clobber to overwrite existing files." )
                return(4)
            if args.clobber and os.path.exists(pdb_file_name):
                print( "File %s already exists -- empty file now" % pdb_file_name )
                open(pdb_file_name, 'w').close()
        for i, molfile in enumerate(molfiles):
            # if the molecule's atom order doent matched with that of the first molfile, skip
            if not compare_molfiles(m, molfile):
                print("conformer %d is not the same molecule with conformer 1, skip" % i)
                break
            #writing out into seperate pdb files
            if not args.all_in_one_pdb:
                pdb_file_name = "%s_%04i.pdb" % (args.pdb, i+1)
                pdb_file = open(pdb_file_name, 'w')
                if not args.clobber and os.path.exists(pdb_file_name):
                    print( "File %s already exists -- skip!" % pdb_file )
                    print( "Use --clobber to overwrite existing files." )
                    break
            else:
                # write into the only 1 pdb file
                pdb_file_name = "%s_rotamer.pdb" % (args.pdb)
                pdb_file = open(pdb_file_name, 'a')

            # m is used for names, molfile is used for XYZ
            write_ligand_pdb(pdb_file, m, molfile, args.name, ctr)
            print("Wrote PDB file %s" % pdb_file_name)

        # Append PDB_ROTAMERS file to params
        if args.use_pdb_rotamers == True and args.use_parent_rotamers is None:
            print("Using PDB rotamers")
            if num_frags == 1: param_file = "%s.params" % args.name
            else: param_file = "%s%i.params" % (args.pdb, i+1)
            write_pdb_rotamers( param_file, pdb_file_name)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

    # Vestigal code for validating automatic atom typing:
    #m = read_mdl_molfile(sys.argv[1])
    #add_fields_to_atoms(m.atoms)
    #assign_rosetta_types(m.atoms)
    #for i, a in enumerate(m.atoms):
    #    err_flag = ""
    #    if a.name.strip() != a.ros_type.strip():
    #        if a.name == "CH1" and a.ros_type.startswith("CH"): pass
    #        else: err_flag = "********" #raise ValueError("typing mismatch!")
    #    print( "%3i %4s --> %4s %s" % (i+1, a.name, a.ros_type, err_flag) )

