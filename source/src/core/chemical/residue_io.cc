// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/residue_io.cc
/// @brief  helper methods for ResidueType input/output
/// @author Phil Bradley

// Unit header
#include <core/chemical/residue_io.hh>

// Package headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/residue_support.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/util.hh>
#include <core/chemical/Bond.hh>
#include <core/chemical/Elements.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecificationFactory.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/DunbrackRotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/PDBRotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/PeptoidRotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/NCAARotamerLibrarySpecification.hh>

// Project headers
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/kinematics/Stub.hh>

// Basic headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/constants.hh>
//#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
//#include <utility/Bound.hh>
#include <utility/vector1.hh>
//#include <utility/keys/AutoKey.hh>
//#include <utility/keys/SmallKeyVector.hh>
#include <utility/io/izstream.hh>

// External headers
#include <boost/foreach.hpp>
#include <boost/graph/graphviz.hpp>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray2D.hh>

// C++ headers

#include <sstream>

namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS


namespace core {
namespace chemical {

static THREAD_LOCAL basic::Tracer tr( "core.chemical" );


///////////////////////////////////////////////////////////////////////////////
/// @brief helper fxn
id::AtomID
atom_id_from_icoor_line(
	std::string const & name,
	ResidueType const & rsd
)
{
	using id::AtomID;
	ICoorAtomID id( name, rsd );

	switch ( id.type() ) {
	case ICoorAtomID::INTERNAL :
		return AtomID( id.atomno(), 1 );
	case ICoorAtomID::CONNECT :
		return AtomID( id.atomno(), 2 );
	case ICoorAtomID::POLYMER_LOWER :
		return AtomID( 1, 3 );
	case ICoorAtomID::POLYMER_UPPER :
		return AtomID( 2, 3 );
	default :
		utility_exit_with_message( "unrecognized stub atom id type!" );
	}
	return id::BOGUS_ATOM_ID;
}


// If polymer, determine a list of main chain atoms by shortest path from LOWER to UPPER.
AtomIndices
define_mainchain_atoms( ResidueTypeOP rsd )
{
	AtomIndices mainchain;

	if ( rsd->is_polymer() ) {
		// Test that this is really a polymer residue.
		if ( rsd->upper_connect_id() && rsd->upper_connect_atom() &&
				rsd->lower_connect_id() && rsd->lower_connect_atom() ) {
			Size upper_connect( rsd->upper_connect_atom() ), lower_connect( rsd->lower_connect_atom() );
			// Default main chain: defined by shortest path from LOWER to UPPER.
			// IMO, everyone should really explicitly define the main chain from the topology file.  ~Labonte
			FArray2D_int D( get_residue_path_distances( *rsd ) );
			uint atom( lower_connect );
			while ( atom != upper_connect ) {
				mainchain.push_back( atom );
				AtomIndices const & nbrs( rsd->nbrs( atom ) );
				int min_d( D( atom, upper_connect ) );
				uint next_atom( atom );

				for ( uint i=1; i<= nbrs.size(); ++i ) {
					uint const nbr( nbrs[i] );
					if ( D( nbr, upper_connect ) < min_d ) {
						min_d = D( nbr, upper_connect );
						next_atom = nbr;
					}
				}
				debug_assert( next_atom != atom );
				atom = next_atom;
			}
			mainchain.push_back( upper_connect );
		} else {
			tr.Warning << "WARNING: Residue " << rsd->name() << " claims it's a polymer, " <<
				"but it doesn't have the appropriate UPPER and LOWER connection points specified.  " <<
				"Set MAINCHAIN_ATOMS in the topology file to remove this warning." << std::endl;
		}
	}
	return mainchain;
}


ResidueTypeOP
read_topology_file(
	std::string const & filename,
	chemical::AtomTypeSetCAP atom_types,
	chemical::ElementSetCAP elements,
	chemical::MMAtomTypeSetCAP mm_atom_types,
	chemical::orbitals::OrbitalTypeSetCAP orbital_atom_types,
	// chemical::CSDAtomTypeSetCAP csd_atom_types kwk commenting out csd_atom_types until I have a chance to fully implement them.
	chemical::ResidueTypeSetCAP rsd_type_set
)
{
	std::string full_filename = filename;
	if ( ! utility::file::file_exists( full_filename ) ) {
		full_filename =  basic::database::full_name( "chemical/residue_type_sets/fa_standard/residue_types/"+filename );
		if ( ! utility::file::file_exists( full_filename ) ) {
			utility_exit_with_message("Cannot find file '"+filename+" or "+full_filename+"'");
		}
	}
	utility::io::izstream data( full_filename.c_str() );
	if ( !data.good() ) {
		utility_exit_with_message("Cannot open file '"+full_filename+"'");
	}
	return read_topology_file(data, atom_types, elements, mm_atom_types, orbital_atom_types, rsd_type_set);
}

///////////////////////////////////////////////////////////////////////////////
/// @details Construct a ResidueType from a file. Example files are currently in
///  main/database/chemical/residue_type_sets/fa_standard/residue_types/l-caa/ directory
///  These files contain information about each basic ResidueType which can be
///  patched to created various variant types.
///
/// The topology file (.params file) is formatted as follows:
///
/// The file may contain any number of lines.  Blank lines and lines beginning with "#"
/// are ignored. Each non-ignored line must begin with a string, a "tag", which
/// describes a piece of data for the ResidueType.  The tags may be given in any order,
/// though they will be processed so that ATOM tag lines are read first.
///
/// Valid tags are:
/// AA:
/// Gives the element of the AA enumeration (src/core/chemical/AA.hh) that
/// is appropriate for this residue type.  This information is used by
/// the knowledge-based potentials which already encode information
/// specifically for proteins or nucleic acids, and is also used by
/// the RotamerLibrary to retrieve the appropriate rotamer library.
/// Provide "aa_unk" here for "unknown" if not dealing with a canonical
/// amino or nucleic acid.  E.g., "AA SER" from SER.params
///
/// ACT_COORD_ATOMS:
/// Lists the atoms that define the "action coordinate" which is used
/// by the fa_pair potential followed by the "END" token. E.g.,
/// "ACT_COORD_ATOMS OG END" from SER.params.
///
/// ADDUCT:
/// Defines an adduct as part of this residue type giving: a) the name, b)
/// the adduct atom name, c) the adduct atom type, d) the adduct mm type,
/// e) the adduct partial charge, and the f) the distance, g) improper bond
/// angle, and h) dihedral that describe how to build the adduct from the
/// i) parent, i) grandparent, and j) great grandparent atoms. E.g.,
/// "ADDUCT  DNA_MAJOR_GROOVE_WATER  WN6 HOH H 0.0   -6.000000  44.000000     2.990000   N6    C6    C5"
/// from ADE.params.
///
/// ATOM:
/// Declare a new atom by name and list several of its properties.
/// This line is column formatted.  The atom's name must be in columns
/// 6-9 so that "ATOM CA  ..." declares a different atom from
/// "ATOM  CA ...".  This is for PDB formatting.  All atom names
/// must be distinct, and all atom names must be distinct when ignoring
/// whitespace ("CA  " and " CA " could not coexist). After the atom name
/// is read, the rest of the line is simply whitespace delimited.
/// Next, the (Rosetta) atom type name is given (which must have
/// been defined in the input AtomTypeSet), and then the mm atom type
/// name is given (which must have been defined in the input
/// MMAtomTypeSet).  Finally, the charge for this atom is given, either
/// as the next input or (ignoring the next input) the one after,
/// if the "parse_charge" flag is on (whatever that is).
/// E.g., "ATOM  CB  CH3  CT3  -0.27   0.000" from ALA.params.
///
/// ATOM_ALIAS:
/// Add alternative name(s) for a given atom, for example to be used when
/// loading a PDB file. This line is column formatted. The canonical
/// Rosetta atom name should be in columns 12-15, with the alternative
/// names coming in columns 17-20, 22-25, etc: "ATOM_ALIAS RRRR 1111 2222 3333 ..."
/// As with the ATOM line, whitespace matters. and aliases must be unique with
/// respect to each other and with canonical names, even when whitespace is ignored
///
/// BACKBONE_AA:
/// Sets the "backbone_aa" for a particular residue, which can be used
/// to template the backbone scoring (rama and p_aa_pp terms).  For example,
/// "BACKBONE_AA ILE" in the non-canonical 4,5-dihydroxyisoleucine params file
/// tells Rosetta to use isoleucine's ramachandran map and p_aa_pp scoring for
/// this noncanonical.
///
/// BOND:
/// Declares a bond between two atoms giving their names. This line is
/// whitespace delimited.  E.g., "BOND  N    CA" from ALA.params.
///
/// BOND_TYPE:
/// Declares a bond between two atoms, giving their names, and also
/// describing the chemical nature of the bond. (The BOND_TYPE line takes the place of
/// a BOND line - do not specify both.) Standard SDF-style numeric
/// descriptors are accepted. (1, 2, 3 for single, double, triple), along with
/// text version SINGLE, DOUBLE, TRIPLE, UNK (unknown), PSEUDO (pseudo bond), ORBITAL,
/// ARO, AMIDE, CARBOXY (for delocalized carboxylate) and DELOCALIZED. Currently UNK/PSEUDO
/// are treated identically, as are ARO/AMIDE/CARBOXY/DELOCALIZED.
/// See convert_to_BondName() in src/core/chemical/Bond.cc for details on parsing.
///
/// CHARGE:
/// Declares a charge for a particular atom.
/// Format CHARGE atom type value
/// Currently valid types are FORMAL. (Partial charges are handled above.)
/// E.g. "CHARGE OG2 FORMAL -1"
///
/// CHI:
/// A misnomer for non-amino acids, declares a side-chain dihedral, its index, and the four atoms that define it.
/// E.g., "CHI 2  CA   CB   CG   CD1" from PHE.params.
///
/// CHI_ROTAMERS:
/// Lists the chi mean/standard-deviation pairs that define how to build
/// rotamer samples.  This is useful for residue types which do not
/// come with their own rotamer libraries.  E.g., "CHI_ROTAMERS 2 180 15"
/// from carbohydrates/to5-beta-D-Psip.params.
///
/// CONNECT:
/// Declares that an inter-residue chemical bond exists from a given
/// atom.  E.g. "CONNECT SG" from CYD.params.  NOTE: Connection order
/// is assumed to be preserved between residue types: connection #2 on
/// one residue type is assumed to be "the same" as connection #2
/// on another residue type, if the two residue types are going
/// to be exchanged in the packer (as ALA might be swapped out
/// for ARG).  CONNECT tags are processed in the order they are
/// listed in the input file.  For polymeric residue types
/// (e.g., proteins, DNA, RNA, saccharides) "LOWER_CONNECT" and "UPPER_CONNECT"
/// should be listed before any additional CONNECT records.  CONNECTs are
/// assigned an index beginning after the LOWER_CONNECT and UPPER_CONNECT, if
/// present.  That is, if a topology file lists both a LOWER_CONNECT and an UPPER_CONNECT,
/// the 1st CONNECT will we given the index 3.
///
/// CUT_BOND:
/// Declares a previously-declared bond to be off-limits to the
/// basic atom-tree construction logic (user-defined atom trees
/// can be created which defy this declaration, if desired).
/// This is useful in cases where the chemical graph contains
/// cycles.  E.g. "CUT_BOND O4' C1'" from URA.params.
///
/// FIRST_SIDECHAIN_ATOM:
/// Gives the name of the first side-chain atom.  All heavy atoms that were
/// declared before the first side-chain atom in the topology file
/// are considered backbone atoms (but not necessarily main-chain atoms).
/// All heavy atoms after the first side-chain atom are considered side-chain atoms.
/// Hydrogen atoms are either side-chain or backbone depending on the heavy atom to
/// which they are bound.  E.g., "FIRST_SIDECHAIN_ATOM CB" from SER.params.
///
/// IO_STRING:
/// Gives the three-letter and one-letter codes that are used to read and
/// write this residue type from and to PDB files, and to FASTA files.
/// This tag is column formatted.  Columns 11-13 are for the three-letter
/// code. Column 15 is for the 1-letter code.  E.g., "IO_STRING Glc Z".
///
/// INTERCHANGEABILITY_GROUP:
/// Gives the name for the group of ResidueType objects that are functionally
/// Interchangeable (but which may have different variant types).  This
/// information is used by the packer to discern what ResidueType to put
/// at a particular position.  If this tag is not given in the topology file,
/// the three-letter code given by the IO_STRING tag is used instead.
///
/// ICOOR_INTERNAL:
/// Describes the geometry of the residue type from internal coordinates
/// giving a) the atom, b) the torsion, phi, in degrees c) the improper bond angle
/// that is (180-bond angle) in degrees, theta, d) the bond length, d, in Angstroms
/// e) the parent atom, f) the grand-parent, and g) the great-grandparent.
/// The first three atoms in the file have a peculiar structure where:
/// 1) at1 lists itself as its own parent, at2 as its grand parent, and at3 as its great-grandparent,
/// 2) at2 lists at1 as its parent, itself as its grand parent, and at3 as its great-grandparent, and
/// 3) at3 list at2 as its parent, at1 as its grand parent, and itself as its great-grandparent.
/// The atoms "LOWER" and "UPPER" are given for polymeric residues to describe the
/// ideal coordinate of the previous and next residues. For non-polymeric inter-residue
/// connections, atoms "CONN#" should be given (e.g., CONN3 for the disulfide connection in CYD).
/// The number for an inter-residue connection comes from the order in which the connection is
/// declared in the file, and includes the LOWER_CONNECT and UPPER_CONNECT connections in this
/// count (e.g., for CYD, there is a LOWER_CONNECT, and UPPER_CONNECT, and only a single
/// CONNECT declaration, so the disulfide connection is the third connection). The order in which
/// internal coordinate data for atoms are given, excepting the first three, must define a "tree"
/// in that atom geometry must only be specified in terms of atoms whose geometry has already been
/// specified.  Improper dihedrals may be specified, where the great grandparent is not the parent atom
/// of the grandparent but in these cases, the great grandparent does need to be a
/// child of the grandparent. E.g.,
/// "ICOOR_INTERNAL    CB  -122.000000   69.862976    1.516263   CA    N     C" from SER.params.
///
/// LOWER_CONNECT:
/// For a polymer residue, declares which atom forms the "lower" inter-residue
/// connection (chemical bond), i.e., the bond to residue i-1.  E.g.,
/// "LOWER_CONNECT N" from SER.params.
///
/// LOWEST_RING_CONFORMER:
/// For a cyclic residue, declares which ideal ring conformation is most stable by IUPAC name.
/// There is no check in place for valid IUPAC conformer names; if the name is not found in the database, no lowest
/// conformer will be set.
/// E.g., "LOWEST_RING_CONFORMER  4C1" from to3-alpha-D-Glcp.params.
///
/// LOW_RING_CONFORMERS:
/// For a cyclic residue, declares which ideal ring conformations are local minima/stable by IUPAC name.
/// (If present, the LOWEST_RING_CONFORMER will automatically included to this subset within RingConformerSet.)
/// There is no check in place for valid IUPAC conformer names; if the name is not found in the database, that conform-
/// er will not be added to the subset.
/// E.g., "LOW_RING_CONFORMERS  O3B B14 3S1 5S1 2SO BO3 1S3 14B 1S5 B25 OS2 1C4" from to3-alpha-D-Glcp.params.
///
/// MAINCHAIN_ATOMS:
/// This is a list of atom names that define the main chain.  The main chain describes the linear connection of atoms
/// from the lower-terminus to the upper-terminus in a residue.  This is NOT synonymous with "backbone atoms".
/// (Backbone atoms are any atoms NOT included in a side chain, as defined by FIRST_SIDECHAIN_ATOM.  See above.)  All
/// main-chain atoms will necessarily be backbone atoms, but not all backbone atoms are main-chain atoms because some
/// residues include rings and/or non-rotatable functional groups.  For example, the carbonyl oxygen of an amino acid
/// residue is a backbone atom but NOT a part of the main chain.
/// If a topology file does not include a MAINCHAIN_ATOMS record, Rosetta will determine the main chain by finding the
/// shortest path from lower terminus to upper terminus, which may be through any CUT_BONDs you have defined!
/// Use of this tag is required for those residue types that ONLY come in LOWER_TERMINUS or UPPER_TERMINUS varieties,
/// such as any residue type that serves exclusively as a "cap" for a larger polymer.
/// E.g., "MAINCHAIN_ATOMS  C1 C2 C3 C4 O4" from an aldohexopyranose topology file.
///
/// METAL_BINDING_ATOMS:
/// For polymer residue types that can bind metals (METALBINDING property), this is a list of the atoms
/// that can form a direct bond to the metal.  For example, in ASP.params, it would read:
/// "METAL_BINDING_ATOMS OD1 OD2"
/// DISULFIDE_ATOM_NAME:
/// For polymer residue types that can form disulfide bonds (SIDECHAIN_THIOL property), this is the atom it does so with.  For example, in CYS.params, it would read:
/// "DISULFIDE_ATOM_NAME SG"
///
/// NAME:
/// Gives the name for this ResidueType.  The name for each ResidueType
/// must be unique within a ResidueTypeSet.  It is not limited to three letters.
/// E.g., "NAME SER" from SER.params.
///
/// NBR_ATOM:
/// Declares the name of the atom which will be used to define the
/// center of the residue for neighbor calculations.  The coordinate
/// of this atom is used along side the radius given in in the
/// NBR_RADIUS tag.  This atom should be chosen to not move
/// during packing so that neighbor relationships can be discerned
/// for an entire set of rotamers from the backbone coordinates from
/// the existing residue.  E.g., "NBR_ATOM CB" from SER.params.
///
/// NBR_RADIUS:
/// Declares a radius that defines a sphere which, when centered
/// on the NBR_ATOM, is guaranteed to contain all heavy atoms
/// under all possible dihedral angle assignments (but where
/// bond angles and lengths are considered ideal).  This is
/// used to determine which residues are close enough that they
/// might interact.  Only the interactions of those such residues
/// are ever evaluated by the scoring machinery.  E.g.,
/// "NBR_RADIUS 3.4473" from SER.params.
///
/// NCAA_ROTLIB_PATH:
/// Gives the path to the rotamer library file to use for a non-canonical
/// amino acid.  E.g., "NCAA_ROTLIB_PATH E35.rotlib" from
/// d-ncaa/d-5-fluoro-tryptophan.params
/// See the ROTAMERS tag for a more general way of specifying how to
/// build rotamers.
///
/// NCAA_ROTLIB_NUM_ROTAMER_BINS:
/// Lists the number of rotamers and the number of bins for each rotamer.
/// E.g., "NCAA_ROTLIB_NUM_ROTAMER_BINS 2 3 2" from
/// d-ncaa/d-5-fluoro-tryptophan.params
/// See the ROTAMERS tag for a more general way of specifying how to
/// build rotamers.
///
/// NCAA_SEMIROTAMERIC:
/// Indicates if a NCAA is semirotameric (this is hard coded for canonicals)
///
/// NRCHI_START_ANGLE:
/// The lower bound for non rotameric chi sampling (sometimes you don't want
/// 0 or 180 to be minimum so you get both sides of these critical values;
/// take this from the input file you used with MakeRotLib
///
/// NRCHI_SYMMETRIC:
/// Indicates if the nonrotameric chi is symmetric (PHE); if absent it is
/// assumed not (ASN)
///
/// NU:
/// Declares an internal ring dihedral, its index, and the four atoms that define it.
/// E.g., "NU 2  C1   C2   C3   C4 ".
///
/// NUMERIC_PROPERTY:
/// Stores an arbitrary float value that goes along with an arbitrary
/// string key in a ResidueType.  No examples can be currently
/// found in the database (10/13).  E.g., "NUMERIC_PROPERTY twelve 12.0"
/// would be a way to store the number "12.0" with the key "twelve".
///
/// ORIENT_ATOM:
/// Describes how to orient rotamers onto an existing residue either
/// by orienting onto the NBR_ATOM atom, or by using the more
/// complicated (default) logic in ResidueType::select_orient_atoms.
/// There are two options here: "ORIENT_ATOM NBR" (orient onto NBR_ATOM) and
/// "ORIENT_ATOM DEFAULT". If this tag is not given in the topology
/// file, then the default behavior is used.  E.g., "SET_ORIENT_ATOM NBR" from
/// SC_Fragment.txt
///
/// PDB_ROTAMERS:
/// Gives the file name that describes entire-residue rotamers which
/// can be used in the packer to consider alternate conformations for
/// a residue.  This is commonly used for small molecules.
/// See the ROTAMERS tag for a more general way of specifying how to
/// build rotamers.
///
/// PROPERTIES:
/// Adds a given set of property strings to a residue type.  E.g.,
/// "PROPERTIES PROTEIN AROMATIC SC_ORBITALS" from TYR.params.
///
/// PROTON_CHI:
/// Declares a previously-declared chi angle to be a "proton chi"
/// and describe how this dihedral should be sampled by code
/// that takes discrete sampling of side-chain conformations.
/// The structure of these samples is as follows:
/// First the word "SAMPLES" is given.  Next the number of samples
/// should be given.  Next a list of dihedral angles, in degrees, should
/// be given.  Next, the word "EXTRA" is given.  Next the number
/// of extra samples is given.  Finally, a list of perturbations
/// angles, in degrees, is given.  In cases where extra rotamers
/// are requested (e.g., with the -ex2 flag), then the listed samples
/// are are perturbed +/- the listed perturbations.  E.g.,
/// "PROTON_CHI 2 SAMPLES 18 0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 EXTRA 0"
/// from SER.params.
///
/// REMAP_PDB_ATOM_NAMES:
/// When reading in a PDB, attempt to match the input atoms for this residue
/// based on elements and estimated connectivity, rather than atom names.
/// (Connectivity by CONECT lines is ignored.) This only applies to input poses,
/// and not to individually loaded residues, e.g. from the PDB_ROTAMERS line.
///
/// ROTAMER_AA:
/// Sets the "rotamer_aa" for a particular residue, which can be used
/// to describe to the RotamerLibrary what amino acid to mimic for the
/// sake of building rotamers.  E.g., "ROTAMER_AA SER"
/// See the ROTAMERS tag for a more general way of specifying how to
/// build rotamers.
///
/// ROTAMERS:
/// Sets the rotamer generation to a given type. The line should be formatted as
/// "ROTAMERS <TAG> <optional data>", where <TAG> specifies what sort of rotamer
/// library to build for this residue, and <optional data> is type-specific data
/// describing how to build.
/// E.g. "ROTAMERS DUNBRACK SER" will cause Serine-like Dunbrack rotamers
/// to be used for this residue type, and "ROTAMERS PDB /path/to/ligand.conf.pdb"
/// will cause the rotamer library to be read from the PDB file "/path/to/ligand.conf.pdb"
/// Currently the following types are known (more may be added):
/// * BASIC - simple rotamer library with no parameters (Just the input residue)
/// * DUNBRACK - Dunbrack rotamer library for the given aa (e.g. "SER")
/// * CENROT - Centroid rotamers for the given aa (e.g. "SER")
/// * NCAA - Non-canonical amino acid libraries (<filename> <bins for each rotamer>)
/// * PEPTOID - Peptoid libraries (<filename> <bins for each rotamer>)
/// * PDB - Rotamer library loaded from PDB (<filename>) - Note that unlike PDB_ROTAMERS,
///      the filename is not relative to the params file location.
///
/// STRING_PROPERTY:
/// Stores an arbitrary string value with a given string key.
/// No example can be currently found in the database (10/13).
/// A valid case would be "STRING_PROPERTY count twelve" which
/// could store the string "twelve" for the key "count".
///
/// TYPE:
/// States whether this is a polymeric or ligand residue type.
/// E.g., "TYPE POLYMER" or "TYPE LIGAND" which adds either
/// "POLYMER" or "LIGAND" properties to this residue type.
///
/// UPPER_CONNECT:
/// For a polymer residue, declares which atom forms the "upper" inter-residue
/// connection (chemical bond), i.e., the bond to residue i+1.  E.g.,
/// "UPPER_CONNECT C" from SER.params.
///
/// VARIANT:
/// Declares this residue type to have a particular variant type.
/// Variant types are used by the packer to determine which
/// ResidueTypes are compatible with a given starting residue.
/// Variants are similar to properties, except that the packer
/// does not restrict itself to residue types that have the
/// same set of properties.  Variant information is also
/// used by the residue-type-patching system to avoid applying
/// patches to certain residue types.  E.g., "VARIANT DISULFIDE".
/// from CYD.params.
///
/// VIRTUAL_SHADOW:
/// Declares the first atom as a shadower of the second atom, implying
/// that the atoms ought to be restrained to lie directly on top of each
/// other. E.g. "VIRTUAL_SHADOW NV N" from PRO.params.  Currently, the
/// cart_bonded and ring_close energy terms are the only energy terms
/// that enforce this.
ResidueTypeOP
read_topology_file(
	utility::io::izstream & data,
	chemical::AtomTypeSetCAP atom_types,
	chemical::ElementSetCAP elements,
	chemical::MMAtomTypeSetCAP mm_atom_types,
	chemical::orbitals::OrbitalTypeSetCAP orbital_atom_types,
	//chemical::CSDAtomTypeSetCAP csd_atom_types kwk commenting out until they have been fully implemented
	chemical::ResidueTypeSetCAP rsd_type_set_ap )
{

	using id::AtomID;
	using id::DOF_ID;
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	using namespace basic;

	std::string filename = data.filename();

	// read the file
	std::string line;
	utility::vector1< std::string > lines;

	std::string myname;
	while ( getline( data, line ) ) {
		//if ( line.size() < 1 || line[0] == '#' ) continue;
		if ( line.size() < 1 ) continue;
		std::string::size_type pound = line.find('#', 0);
		if ( pound == std::string::npos ) {
			lines.push_back( line );
		} else {
			std::string no_comment_line= line.substr(0, pound);
			lines.push_back(no_comment_line);
		}
		if ( line.size() > 5 && line.substr(0,5)=="NAME " ) {
			std::istringstream l( line );
			std::string tag;
			l >> tag >> myname;
			if ( l.fail() ) myname.clear();
		}
	}
	data.close();

	std::map< std::string, std::string > atom_type_reassignments;
	std::map< std::string, Real > atomic_charge_reassignments;
	if ( rsd_type_set_ap.lock() ) {
		setup_atom_type_reassignments_from_commandline( myname, rsd_type_set_ap.lock()->name(), atom_type_reassignments );
		setup_atomic_charge_reassignments_from_commandline( myname, rsd_type_set_ap.lock()->name(), atomic_charge_reassignments );
	}

	// Decide what type of Residue to instantiate.
	// would scan through for the TYPE line, to see if polymer or ligand...
	//
	// Residue needs a pointer to the AtomTypeSet object for setting up atom-type dependent data.
	//
	// You may be asking yourself, at least I was, why the hell do we scan this file more than once?
	// The problem is that while reading the params (topology file), we need to assign private member variable
	// data in ResidueType.  We want to make sure that we get the correct number of atoms assigned
	// before we start adding bonds, icoor, etc.  This allows to provide checks to make sure certain
	// things are being assigned correctly, i.e., adding bonds correctly, setting icoor values with correct placement
	// of stub atoms, etc., etc.

	ResidueTypeOP rsd( new ResidueType( atom_types.lock(), elements.lock(), mm_atom_types.lock(), orbital_atom_types.lock() ) ); //kwk commenting out until atom types are fully implemented , csd_atom_types ) );
	if ( ! rsd_type_set_ap.expired() ) {
		rsd->residue_type_set( rsd_type_set_ap );  // Give this rsd_type a backpointer to its set.
	}

	// Add the atoms.
	Size const nlines( lines.size() );
	Size natoms(0);//, norbitals(0);
	for ( Size i=1; i<= nlines; ++i ) {
		std::string line( lines[i] );
		if ( line.size() > 0 ) {
			while ( line.substr(line.size()-1) == " " ) {
				line = line.substr(0, line.size()-1);
				if ( line.size() == 0 ) break;
			}
		}

		std::istringstream l( line );
		std::string tag;

		l >> tag;
		// if ( line.size() < 5 || line.substr(0,5) != "ATOM " ) continue;
		if ( tag != "ATOM" ) continue;

		// the atom name for this atom
		std::string const atom_name( line.substr(5,4) );
		l >> tag; // std::string const atom_name( tag );

		// the atom type name -- must match one of the atom types
		// for which force-field parameters exist
		// std::string const atom_type_name( line.substr(10,4) );
		l >> tag; std::string atom_type_name( tag ); // non-const in case of atomtype reassignment

		// read in the Molecular mechanics atom type
		// std::string const mm_atom_type_name( line.substr(15,4) );
		l >> tag; std::string const mm_atom_type_name( tag );

		// the atomic charge
		float charge;
		// std::istringstream l( line.substr(20) );
		l >> charge;
		float parse_charge(charge);
		if ( !l.eof() ) {
			l >> parse_charge;
		}

		utility::vector1< std::string > props;
		std::string prop;
		while ( l >> prop ) {
			props.push_back( prop );
		}

		if ( atom_type_reassignments.find( stripped( atom_name ) ) != atom_type_reassignments.end() ) {
			tr.Trace << "reassigning atom " << atom_name << " atomtype: " << atom_type_name << " --> " <<
				atom_type_reassignments.find( stripped( atom_name ) )->second << ' ' << filename << std::endl;
			atom_type_name = atom_type_reassignments.find( stripped( atom_name ) )->second;
		}

		if ( atomic_charge_reassignments.find( stripped( atom_name ) ) != atomic_charge_reassignments.end() ) {
			tr.Trace << "reassigning atomic charge " << atom_name << " atomtype: " << atom_type_name << " --> " <<
				atomic_charge_reassignments.find( stripped( atom_name ) )->second << ' ' << filename << std::endl;
			// note that we set charge and also parse_charge, so this will over-ride the parse_charge if those are the charges we are using
			charge = parse_charge = atomic_charge_reassignments.find( stripped( atom_name ) )->second;
		}

		if ( ! basic::options::option[ basic::options::OptionKeys::corrections::chemical::parse_charge ]() ) {
			rsd->add_atom( atom_name, atom_type_name, mm_atom_type_name, charge );
		} else {
			rsd->add_atom( atom_name, atom_type_name, mm_atom_type_name, parse_charge );
		}

		for ( Size p = 1; p <= props.size(); ++p ) {
			rsd->atom( atom_name ).set_property( props[ p ], true );
		}

		++natoms;
	}

	// No ATOM lines probably means an invalid file.  Perhaps someone made a mistake with an -extra_res_fa flag.
	// Fail gracefully now, versus a segfault later.
	if ( natoms == 0 ) {
		utility_exit_with_message("Residue topology file '" + filename + "' does not contain valid ATOM records.");
	}

	// Add the bonds; parse the rest of file.
	bool found_AA_record = false;
	AtomIndices mainchain_atoms;

	// Set disulfide atom name to "NONE"
	// So that's the default
	rsd->set_disulfide_atom_name( "NONE" );

	for ( Size i=1; i<= nlines; ++i ) {
		std::string const & line( lines[i] );
		std::istringstream l( line );
		std::string tag,atom1,atom2,atom3,atom4, rotate, /*orbitals_tag, orbital,*/ bond_type;
		std::string pdb_rotamers_filename;
		core::Real value;
		l >> tag;
		if ( l.fail() ) continue;
		if ( tag == "CONNECT" ) {
			l >> atom1;
			l >> rotate; // not used here
			rsd->add_residue_connection( atom1);
		} else if ( tag == "TYPE" ) {
			// will probably handle this differently later on
			l >> tag;
			if ( tag == "POLYMER" ) {
				rsd->add_property( tag );
			} else if ( tag == "LIGAND" ) {
				rsd->add_property( tag );
			}
		} else if ( tag == "METAL_BINDING_ATOMS" ) {
			l >> atom1;
			while ( !l.fail() ) {
				rsd->add_metalbinding_atom( atom1 );
				l >> atom1;
			}
		} else if ( tag == "DISULFIDE_ATOM_NAME" ) {
			l >> tag;
			rsd->set_disulfide_atom_name( tag );
		} else if ( tag == "ATOM_ALIAS" ) {
			if ( line.size() < 20 ) {
				utility_exit_with_message("ATOM_ALIAS line too short in " + filename + ":\n" + line );
			}
			atom1 = line.substr( 11, 4 ); // Rosetta atom
			if ( ! rsd->has( atom1 ) ) {
				utility_exit_with_message( "ATOM_ALIAS line in " + filename + " attempts to add atom " + atom1 + " which is not a member of " + rsd->name() );
			}
			core::Size pos(16);
			while ( line.size() >= pos+4 ) {
				atom2 = line.substr(pos, 4);
				rsd->add_atom_alias( atom1, atom2 );
				pos += 5;
			}
		} else if ( tag == "BOND" ) {
			l >> atom1 >> atom2;
			rsd->add_bond( atom1, atom2 );
		} else if ( tag == "BOND_TYPE" ) {
			l >> atom1 >> atom2 >> bond_type;
			rsd->add_bond(atom1, atom2, convert_to_BondName(bond_type));
		} else if ( tag == "CHARGE" ) {
			l >> atom1;
			// We should allow multiple charges on one line, but since we now just have the one, hold off.
			l >> tag;
			if ( tag == "FORMAL" ) {
				l >> value;
				rsd->atom( atom1 ).formal_charge( int(value) );
			} else {
				utility_exit_with_message("ERROR: Invalid charge type '"+tag+"' in topology file, " + filename );
			}
		} else if ( tag == "CUT_BOND" ) {
			l >> atom1 >> atom2;
			rsd->add_cut_bond( atom1, atom2 );
		} else if ( tag == "CHI" ) {
			Size chino;
			l >> chino >> atom1 >> atom2 >> atom3 >> atom4;
			rsd->add_chi( chino, atom1, atom2, atom3, atom4 );
		} else if ( tag == "NU" ) {
			uint nu_num;
			l >> nu_num >> atom1 >> atom2 >> atom3 >> atom4;
			rsd->add_nu(nu_num, atom1, atom2, atom3, atom4);
		} else if ( tag == "ADD_RING" ) {
			uint ring_num;
			l >> ring_num;
			utility::vector1< std::string > ring_atoms;
			l >> atom1;
			while ( !l.fail() ) {
				ring_atoms.push_back( atom1 );
				l >> atom1;
			}
			rsd->add_ring( ring_num, ring_atoms );
		} else if ( tag == "LOWEST_RING_CONFORMER" ) {
			uint ring_num;
			std::string conformer;
			l >> ring_num >> conformer;
			rsd->set_lowest_energy_ring_conformer( ring_num, conformer );
		} else if ( tag == "LOW_RING_CONFORMERS" ) {
			uint ring_num;
			utility::vector1< std::string > conformers;
			l >> ring_num;
			l >> tag;
			while ( ! l.fail() ) {
				conformers.push_back( tag );
				l >> tag;
			}
			rsd->set_low_energy_ring_conformers( ring_num, conformers );
		} else if ( tag == "PROTON_CHI" ) {
			Size chino, nsamples(0), nextra_samples(0);
			std::string dummy;
			l >> chino;
			l >> dummy; // should be "SAMPLES"
			l >> nsamples;
			utility::vector1< Real > samples( nsamples );
			for ( Size ii = 1; ii <= nsamples; ++ii ) {
				l >> samples[ ii ];
			}
			l >> dummy; // should be "EXTRA"
			l >> nextra_samples;
			utility::vector1< Real > extra_samples( nextra_samples );
			for ( Size ii = 1; ii <= nextra_samples; ++ii ) {
				l >> extra_samples[ ii ];
			}
			if( ! l ) {
				tr << "BAD PROTON_CHI line: " << line << std::endl;
				utility_exit_with_message("Malformed PROTON_CHI line in params file " + filename);
			}
			if ( basic::options::option[ basic::options::OptionKeys::corrections::chemical::expand_st_chi2sampling ]
					&& (rsd->aa() == aa_ser || rsd->aa() == aa_thr )
					&& rsd_type_set_ap.lock() && rsd_type_set_ap.lock()->name() == FA_STANDARD ) {
				// ugly, temporary hack. change the sampling for serine and threonine chi2 sampling
				// so that proton chi rotamers are sampled ever 20 degrees
				tr << "Expanding chi2 sampling for amino acid " << rsd->aa() << std::endl;
				utility::vector1< Real > st_expanded_samples( 18, 0 );
				for ( Size ii = 1; ii <= 18; ++ii ) {
					st_expanded_samples[ ii ] = (ii-1) * 20;
				}
				samples = st_expanded_samples;
				extra_samples.resize(0);
			}
			rsd->set_proton_chi( chino, samples, extra_samples );

		} else if ( tag == "NBR_ATOM" ) {
			l >> atom1;
			rsd->nbr_atom( atom1 );

		} else if ( tag == "NBR_RADIUS" ) {
			Real radius;
			l >> radius;
			rsd->nbr_radius( radius );
		} else if ( tag == "ORIENT_ATOM" ) {
			l >> tag;
			if ( tag == "NBR" ) {
				rsd->force_nbr_atom_orient(true);
			} else if ( tag == "DEFAULT" ) {
				rsd->force_nbr_atom_orient(false);
			} else {
				utility_exit_with_message("Unknown ORIENT_ATOM mode: " + tag + " in " + filename );
			}
		} else if ( tag == "PROPERTIES" ) {
			l >> tag;
			while ( !l.fail() ) {
				rsd->add_property( tag );
				l >> tag;
			}
		} else if ( tag == "NUMERIC_PROPERTY" ) {
			core::Real value = 0.0;
			l >> tag >> value;
			rsd->add_numeric_property(tag,value);

		} else if ( tag == "STRING_PROPERTY" ) {
			std::string value;
			l >> tag >> value;
			rsd->add_string_property(tag,value);

		} else if ( tag == "VARIANT" ) {
			l >> tag;
			while ( !l.fail() ) {
				rsd->add_variant_type( tag );
				l >> tag;
			}
		} else if ( tag == "MAINCHAIN_ATOMS" ) {
			// Note: Main-chain atoms describe the linear connection of atoms from the lower-terminus to the upper-
			// terminus in a residue.  This is NOT synonymous with "backbone atoms".  (Backbone atoms are any atoms NOT
			// included in a side chain, as defined by FIRST_SIDECHAIN_ATOM.  See below.)  All main-chain atoms will
			// necessarily be backbone atoms, but not all backbone atoms are main-chain atoms because some residues
			// include rings and/or non-rotatable functional groups.  For example, the carbonyl oxygen of an amino acid
			// residue is a backbone atom but NOT a part of the main chain.
			// If a topology file does not include a MAINCHAIN_ATOMS record, Rosetta will determine the main chain by
			// finding the shortest path from lower terminus to upper terminus.
			l >> tag;
			while ( !l.fail() ) {
				mainchain_atoms.push_back( rsd->atom_index( tag ) );
				l >> tag;
			}
		} else if ( tag == "FIRST_SIDECHAIN_ATOM" ) {
			// Note: Atoms are side-chain by default.
			// In this case, the opposite of "side-chain" is "backbone"; the "main chain" will be a subset of the
			// backbone atoms, because some residues may contain rings or non-rotatable functional groups, such as a
			// carbonyl.  The main chain is defined separately.  (See above.)
			l >> tag;
			if ( tag == "NONE" ) {
				// Set all atoms to backbone.
				for ( Size j=1; j<= rsd->natoms(); ++j ) {
					rsd->set_backbone_heavyatom( rsd->atom_name(j) );
				}
			} else if ( rsd->has( tag ) ) {
				for ( Size j=1; j< rsd->atom_index( tag ); ++j ) {
					rsd->set_backbone_heavyatom( rsd->atom_name(j) );
				}
			}

		} else if ( tag == "IO_STRING" ) {
			debug_assert( line.size() >= 15 );
			std::string const three_letter_code( line.substr(10,3) ),
				one_letter_code( line.substr(14,1) );
			rsd->name3( three_letter_code );
			rsd->name1( one_letter_code[0] );
			// Default behavior for interchangeability_group is to take name3
			if ( rsd->interchangeability_group() == "" ) {
				rsd->interchangeability_group( three_letter_code );
			}
		} else if ( tag == "INTERCHANGEABILITY_GROUP" ) {
			// The INTERCHANGEABILITY_GROUP tag is not required; the name3 from the IO_STRING will be
			// used instead.
			l >> tag;
			rsd->interchangeability_group( tag );
		} else if ( tag == "AA" ) {
			l >> tag;
			rsd->aa( tag );
			found_AA_record = true;

		} else if ( tag == "BACKBONE_AA" ) {
			l >> tag;
			rsd->backbone_aa( tag );

		} else if ( tag == "NAME" ) {
			l >> tag;
			rsd->name( tag ); //The name will have variant types appended to it; it will be the unique identifer for a ResidueType.
			rsd->base_name( tag ); //The base name stays the same once set.  It's common to all ResidueTypes that share a base type but differ in their VariantTypes.
		} else if ( tag == "CHI_ROTAMERS" ) {
			Size chino;
			Real mean, sdev;
			l >> chino;
			l >> mean >> sdev;
			while ( !l.fail() ) {
				rsd->add_chi_rotamer( chino, mean, sdev );
				l >> mean >> sdev;
			}
		} else if ( tag == "REMAP_PDB_ATOM_NAMES" ) {
			rsd->remap_pdb_atom_names( true );
		} else if ( tag == "ACT_COORD_ATOMS" ) {
			while ( l ) {
				l >> atom1;
				if ( atom1 == "END" ) break;
				rsd->add_actcoord_atom( atom1 );
			}
		} else if ( tag == "LOWER_CONNECT" ) {
			l >> atom1;
			rsd->set_lower_connect_atom( atom1 );
		} else if ( tag == "UPPER_CONNECT" ) {
			l >> atom1;
			rsd->set_upper_connect_atom( atom1 );
		} else if ( tag == "ADDUCT" ) {
			std::string adduct_name, adduct_atom_name, adduct_atom_type, adduct_mm_type;
			Real adduct_q, adduct_d, adduct_theta, adduct_phi;
			l >> adduct_name >> adduct_atom_name;
			l >> adduct_atom_type >> adduct_mm_type;
			l >> adduct_q >> adduct_phi >> adduct_theta >> adduct_d;
			l >> atom1 >> atom2 >> atom3;
			ObjexxFCL::lowercase(adduct_name);
			Adduct new_adduct( adduct_name, adduct_atom_name,
				adduct_atom_type, adduct_mm_type, adduct_q,
				adduct_phi, adduct_theta, adduct_d,
				atom1, atom2, atom3 );
			rsd->add_adduct( new_adduct );
		} else if ( tag == "ROTAMERS" ) {
			using namespace core::chemical::rotamers;
			if ( rsd->rotamer_library_specification() ) {
				tr.Error << "Found existing rotamer specification " << rsd->rotamer_library_specification()->keyname() << " when attempting to set ROTAMERS specification" << std::endl;
				utility_exit_with_message("Cannot have multiple rotamer specifications in params file, " + filename );
			}
			l >> tag;
			if ( ! l ) { utility_exit_with_message("Must provide rotamer library type in ROTAMERS line, from " + filename ); }
			RotamerLibrarySpecificationOP rls( RotamerLibrarySpecificationFactory::get_instance()->get( tag, l ) ); // Create with remainder of line.
			rsd->rotamer_library_specification( rls );
		} else if ( tag == "ROTAMER_AA" ) {
			using namespace core::chemical::rotamers;
			if ( rsd->rotamer_library_specification() ) {
				tr.Error << "Found existing rotamer specification " << rsd->rotamer_library_specification()->keyname() << " when attempting to set ROTAMERS specification" << std::endl;
				utility_exit_with_message("Cannot have multiple rotamer specifications in params file, " + filename );
			}
			tag = DunbrackRotamerLibrarySpecification::library_name();
			RotamerLibrarySpecificationOP rls( RotamerLibrarySpecificationFactory::get_instance()->get( tag, l ) ); // Create with remainder of line (aa)
			rsd->rotamer_library_specification( rls );
			// TODO: Expand for d-aa and nucleic ?
		} else if ( tag == "PDB_ROTAMERS" ) {
			// Assume name of rotamers file has no path info, etc
			// and is in same directory as the residue parameters file.
			using namespace utility::file;
			l >> pdb_rotamers_filename;
			FileName this_file( filename ), rot_file( pdb_rotamers_filename );
			rot_file.vol( this_file.vol() );
			rot_file.path( this_file.path() );

			if ( rsd->rotamer_library_specification() ) {
				tr.Error << "Found existing rotamer specification " << rsd->rotamer_library_specification()->keyname() << " when attempting to set PDB_ROTAMERS parameters." << std::endl;
				utility_exit_with_message("Cannot have multiple rotamer specifications in params file, " + filename );
			}
			rsd->rotamer_library_specification( core::chemical::rotamers::PDBRotamerLibrarySpecificationOP( new core::chemical::rotamers::PDBRotamerLibrarySpecification( rot_file() ) ) );

			tr.Debug << "Setting up conformer library for " << rsd->name() << std::endl;
		} else if ( tag == "NCAA_ROTLIB_PATH" || tag == "NCAA_SEMIROTAMERIC" || tag == "NCAA_ROTLIB_NUM_ROTAMER_BINS" ||
				tag == "NRCHI_SYMMETRIC" || tag == "NRCHI_START_ANGLE" ) {

			using namespace core::chemical::rotamers;
			NCAARotamerLibrarySpecificationOP ncaa_libspec;
			if ( rsd->rotamer_library_specification() ) {
				NCAARotamerLibrarySpecificationCOP old_libspec( utility::pointer::dynamic_pointer_cast< NCAARotamerLibrarySpecification const >( rsd->rotamer_library_specification() ) );
				if ( ! old_libspec ) {
					tr.Error << "Found existing rotamer specification " << rsd->rotamer_library_specification()->keyname();
					tr.Error << " when attempting to set " << tag << " parameter for NCAA rotamer libraries." << std::endl;
					utility_exit_with_message("Cannot have multiple rotamer specifications in params file, " + filename );
				}
				ncaa_libspec = NCAARotamerLibrarySpecificationOP( new NCAARotamerLibrarySpecification( *old_libspec ) );
			} else {
				ncaa_libspec = NCAARotamerLibrarySpecificationOP( new core::chemical::rotamers::NCAARotamerLibrarySpecification );
			}

			if ( tag == "NCAA_ROTLIB_PATH" ) {
				std::string path;
				l >> path;

				ncaa_libspec->ncaa_rotlib_path( path );
			} else if ( tag == "NCAA_SEMIROTAMERIC" ) {
				ncaa_libspec->semirotameric_ncaa_rotlib( true );
			} else if ( tag == "NCAA_ROTLIB_NUM_ROTAMER_BINS" ) {
				Size n_rots(0);
				utility::vector1<Size> n_bins_per_rot;
				l >> n_rots;
				n_bins_per_rot.resize( n_rots );
				for ( Size i = 1; i <= n_rots; ++i ) {
					Size bin_size(0);
					l >> bin_size;
					n_bins_per_rot[i] = bin_size;
				}
				ncaa_libspec->ncaa_rotlib_n_bin_per_rot( n_bins_per_rot );
			} else if ( tag == "NRCHI_SYMMETRIC" ) {
				// this tag present = true
				ncaa_libspec->nrchi_symmetric( true );
			} else if ( tag == "NRCHI_START_ANGLE" ) {
				Real angle(-180);
				l >> angle;
				ncaa_libspec->nrchi_start_angle( angle );
			} else {
				tr.Error << "Did not expect " << tag << " when reading NCAA rotamer info." << std::endl;
				utility_exit_with_message("Logic error in params file, " + filename + " reading.");
			}

			rsd->rotamer_library_specification( ncaa_libspec );

			// End of NCAA library entries.
		} else if ( tag == "PEPTOID_ROTLIB_PATH" || tag == "PEPTOID_ROTLIB_NUM_ROTAMER_BINS" ) {
			using namespace core::chemical::rotamers;
			PeptoidRotamerLibrarySpecificationOP peptoid_libspec;
			if ( rsd->rotamer_library_specification() ) {
				PeptoidRotamerLibrarySpecificationCOP old_libspec( utility::pointer::dynamic_pointer_cast< PeptoidRotamerLibrarySpecification const >( rsd->rotamer_library_specification() ) );
				if ( ! old_libspec ) {
					tr.Error << "Found existing rotamer specification " << rsd->rotamer_library_specification()->keyname();
					tr.Error << " when attempting to set " << tag << " parameter for peptoid rotamer libraries." << std::endl;
					utility_exit_with_message("Cannot have multiple rotamer specifications in params file, " + filename );
				}
				peptoid_libspec = PeptoidRotamerLibrarySpecificationOP( new PeptoidRotamerLibrarySpecification( *old_libspec ) );
			} else {
				peptoid_libspec = PeptoidRotamerLibrarySpecificationOP( new core::chemical::rotamers::PeptoidRotamerLibrarySpecification );
			}

			if ( tag == "PEPTOID_ROTLIB_PATH" ) {
				std::string path;
				l >> path;
				peptoid_libspec->peptoid_rotlib_path( path );
			} else if ( tag == "PEPTOID_ROTLIB_NUM_ROTAMER_BINS" ) {
				Size n_rots(0);
				utility::vector1<Size> n_bins_per_rot;
				l >> n_rots;
				n_bins_per_rot.resize( n_rots );
				for ( Size i = 1; i <= n_rots; ++i ) {
					Size bin_size(0);
					l >> bin_size;
					n_bins_per_rot[i] = bin_size;
				}
				peptoid_libspec->peptoid_rotlib_n_bin_per_rot( n_bins_per_rot );
			} else {
				tr.Error << "Did not expect " << tag << " when reading peptoid rotamer info." << std::endl;
				utility_exit_with_message("Logic error in params file reading, " + filename );
			}

			rsd->rotamer_library_specification( peptoid_libspec );
			// End of peptoid library entries.
		} else if ( tag == "VIRTUAL_SHADOW" ) {
			std::string shadower, shadowee;
			l >> shadower >> shadowee;
			rsd->set_shadowing_atom( shadower, shadowee );
		} else if ( tag == "ATOM" || tag == "ICOOR_INTERNAL" ) {
			; // ATOM lines handled above, ICOOR_INTERNAL lines handled below
		} else {
			tr.Warning << "WARNING: Ignoring line starting with '" << tag << "' when parsing topology file." << std::endl;
		}

	} // i=1,nlines


	if ( !found_AA_record ) {
		tr.Warning << "No AA record found for " << rsd->name()
			<< "; assuming " << name_from_aa( rsd->aa() ) << std::endl;
	}


	// set icoor coordinates, store information about polymer links
	// also sets up base_atom
	{
		std::map< std::string, utility::vector1< std::string > > icoor_reassignments;
		if ( rsd_type_set_ap.lock() ) {
			setup_icoor_reassignments_from_commandline( myname, rsd_type_set_ap.lock()->name(), icoor_reassignments );
		}

		std::map< std::string, Vector > rsd_xyz;  // The coordinates of each atom in the residue.

		for ( Size i=1; i<= nlines; ++i ) {

			std::string const & line( lines[i] );
			std::istringstream l( line );
			std::string tag, child_atom, parent_atom, angle_atom, torsion_atom;

			Real phi, theta, d;
			l >> tag;

			if ( tag != "ICOOR_INTERNAL" ) continue;

			l >> child_atom >> phi >> theta >> d >> parent_atom >> angle_atom >> torsion_atom;

			bool const symm_gly_corrections( basic::options::option[ basic::options::OptionKeys::score::symmetric_gly_tables ].user() ); //Are we symmetrizing the glycine params file?
			if ( icoor_reassignments.find( child_atom ) != icoor_reassignments.end() ) {
				utility::vector1< std::string > const & new_params( icoor_reassignments.find( child_atom )->second );
				phi   = ObjexxFCL::float_of( new_params[1] );
				theta = ObjexxFCL::float_of( new_params[2] );
				d     = ObjexxFCL::float_of( new_params[3] );
				parent_atom  = new_params[4];
				angle_atom   = new_params[5];
				torsion_atom = new_params[6];
				tr.Trace << "reassigning icoor: " << myname << ' ' << child_atom << ' ' <<
					phi << ' ' << theta << ' ' << d  << ' ' <<
					parent_atom << ' ' << angle_atom<< ' ' << torsion_atom << std::endl;
			}

			if ( symm_gly_corrections && rsd->aa() == core::chemical::aa_gly ) { //If the user has used the -symmetric_gly_tables option, we need to symmetrize the glycine params file.
				apply_symm_gly_corrections( child_atom, phi, theta, d, parent_atom, angle_atom, torsion_atom  );
			}

			phi = radians(phi); theta = radians(theta); // in degrees in the file for human readability

			// This code should probably be extracted to a util function
			if ( natoms > 1 ) {
				// build the Cartesian coords for the new atom:
				if ( child_atom == parent_atom ) {
					if ( ! rsd_xyz.empty() ) {
						utility_exit_with_message("Only the first ICOOR atom in a topology file should list itself as its own parent atom, from file " + filename );
					}
					rsd_xyz[ child_atom ] = Vector( 0.0 );

				} else if ( child_atom == angle_atom ) {
					if ( rsd_xyz.size() != 1 ) {
						utility_exit_with_message("Only the second ICOOR atom in a topology file should list itself as its own angle atom, from file " + filename  );
					}
					if ( ! rsd_xyz.count( parent_atom ) ) {
						utility_exit_with_message("In second ICOOR atom in topology file - parent atom not found, from file " + filename );
					}
					rsd_xyz[ child_atom ] = Vector( d, 0.0, 0.0 );
				} else {
					Vector torsion_xyz;
					if ( child_atom == torsion_atom ) {
						if ( rsd_xyz.size() != 2 ) {
							utility_exit_with_message("Only the third ICOOR atom in a topology file should list itself as its own dihedral atom, from file " + filename );
						}
						if ( ! rsd_xyz.count( parent_atom ) || ! rsd_xyz.count( angle_atom ) ) {
							utility_exit_with_message("In third ICOOR atom in topology file - parent and/or angle atom not found, from file " + filename );
						}
						torsion_xyz = Vector( 1.0, 1.0, 0.0 );
					} else {
						if ( ! ( rsd_xyz.count( parent_atom ) && rsd_xyz.count( angle_atom ) &&
								rsd_xyz.count( torsion_atom ) ) ) {
							utility_exit_with_message("In ICOOR atom line in topology file: " + filename +
								"; reference atoms must be specified in earlier line.  Missing " +
								parent_atom + " or " + angle_atom + " or " + torsion_atom);
						}
						torsion_xyz = rsd_xyz[ torsion_atom ];
					}

					kinematics::Stub const stub( rsd_xyz[ parent_atom ], rsd_xyz[ angle_atom ], torsion_xyz );
					rsd_xyz[ child_atom ] = stub.spherical( phi, theta, d );
				}
			}

			// set atom_base
			if ( child_atom != "UPPER" && child_atom != "LOWER" && child_atom.substr(0,4) != "CONN" ) {
				// atom base only valid for genuine atoms of this residue
				if ( child_atom == parent_atom ) {
					// root of the tree
					if ( natoms == 1 ) {
						rsd->set_atom_base( child_atom, child_atom ); // 1st child of root atom
					} else {
						rsd->set_atom_base( child_atom, angle_atom ); // 1st child of root atom
					}
				} else {
					rsd->set_atom_base( child_atom, parent_atom );
				}
			}

			// set icoor
			rsd->set_icoor(child_atom, phi, theta, d, parent_atom, angle_atom, torsion_atom );

		} // loop over file lines looking for ICOOR_INTERNAL lines


		// fill in the rsd-xyz values
		if ( natoms == 1 ) {
			std::string const name( rsd->atom_name(1) );
			rsd->set_ideal_xyz( name, Vector(0.0) );

		} else {
			// now fill in the icoor values -- in principle the rsd itself could be doing this...
			for ( Size i=1; i<= natoms; ++i ) {
				std::string name( rsd->atom_name(i) );
				strip_whitespace( name );
				debug_assert( rsd_xyz.count( name ) );
				rsd->set_ideal_xyz( name, rsd_xyz[ name ] );
				//rsd->set_xyz( rsd->atom_name(i), atom_tree.xyz( id::AtomID(i,1) ) );
				//rsd->atom(i).xyz( atom_tree.xyz( id::AtomID(i,1) ) );
			}
		}


		// If polymer, fill list of main chain atoms, if not already defined by MAINCHAIN_ATOMS.
		// (This must occur after internal coordinates and connections are set.)
		if ( mainchain_atoms.size() == 0 ) {
			mainchain_atoms = define_mainchain_atoms( rsd );
		}
		rsd->set_mainchain_atoms( mainchain_atoms );

		// Okay, now that we have information about the residue's atomic positions
		// we can re-evaluate if there could be an issue with the specification
		// of L_AA and D_AA.
		// Don't worry about GLY. This is properly handled inside the chirality detection
		// but we just don't want to print the extra warning in 99% of poses!
		// AMW: somehow, staple residues 08A and 08B are being counted as protein. WTF?
		// They have no properties set!
		// Not surprising--SOME of the 08A/08B params are set as protein, so depending
		// on your RTS this will happen.
		if ( rsd->is_protein() && !rsd->is_achiral_backbone() && !rsd->is_l_aa() && !rsd->is_d_aa() ) {

			tr.Warning << "Warning: protein residue " << rsd->name3() << " is not explicitly listed"
				<< " as either L or D in its params file." << std::endl;
			tr.Warning << "To avoid seeing this warning in the future, add \"L_AA\", \"D_AA\", "
				<< "or \"ACHIRAL_BACKBONE\" to the \"PROPERTIES\" line of the params file." << std::endl;

			bool is_l_aa = false;
			bool is_d_aa = false;

			detect_ld_chirality_from_polymer_residue( rsd_xyz, rsd->name3(), is_d_aa, is_l_aa );

			tr.Trace << "Detected chirality " << ( is_l_aa ? "L" : ( is_d_aa ? "D" : "ACHIRAL" ) ) << " from params " << std::endl;

			if ( is_l_aa ) rsd->add_property( "L_AA" );
			if ( is_d_aa ) rsd->add_property( "D_AA" );
			debug_assert( !( rsd->is_l_aa() && rsd->is_d_aa() ) );
			if ( !is_l_aa && !is_d_aa ) rsd->add_property( "ACHIRAL_BACKBONE" ); //Automatically set up achirality, too.
			debug_assert( !( rsd->is_l_aa() && rsd->is_achiral_backbone() ) && !( rsd->is_d_aa() && rsd->is_achiral_backbone() ) ); //Double-check that a residue isn't both chiral and achiral.
		}

		// now also need to store the information about the geometry at the links...

	} // scope

	// calculate any remaining derived data
	rsd->finalize();

	return rsd;
}


/// @details function to write out a topology file given a residue type, can be used to
/// debug on the fly generated residue types. Note: not perfect yet, the enums for
/// the connection types are given in numbers instead of names
void
write_topology_file(
	ResidueType const & rsd,
	std::string filename /*= ""*/
)
{
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	if ( ! filename.size() ) {
		filename = rsd.name() + ".params";
	}

	std::ofstream out( filename.c_str() );

	out << "#rosetta residue topology file \n";
	out << "#version ??? \n";
	out << "#This automatically generated file is not really formatted, but should work, excpet that enums for connection types \n# are given in numbers and not strings. \n";

	// first write out all the general tags
	out << "NAME " << rsd.name() << " \n";
	out << "IO_STRING " << rsd.name3() << " " << rsd.name1() << " \n";
	if ( rsd.is_polymer() ) { out << "TYPE POLYMER \n"; }
	else if ( rsd.is_ligand() ) { out << "TYPE LIGAND \n"; }
	else if ( rsd.is_surface() ) { out << "TYPE SURFACE \n"; }
	out << "AA " << rsd.aa() << " \n";

	// then write out the atoms
	for ( Size i=1; i <= rsd.natoms(); ++i ) {

		std::string atom_out = "ATOM " + rsd.atom_name( i ) + " " + rsd.atom_type( i ).name() + "  ";
		atom_out = atom_out + rsd.mm_atom_type(i).name();
		out << atom_out << " " << rsd.atom(i).charge() << " \n";

	} // atom write out

	if ( rsd.is_polymer() ) {
		if ( !rsd.is_lower_terminus() ) { out << "LOWER_CONNECT " << rsd.atom_name( rsd.lower_connect().atomno() ) << " \n"; }
		if ( !rsd.is_upper_terminus() ) { out << "UPPER_CONNECT " << rsd.atom_name( rsd.upper_connect().atomno() ) << " \n";}
	}

	// then all the bonds
	for ( Size i=1; i <= rsd.natoms(); ++i ) {
		BOOST_FOREACH ( Size atom_index, rsd.nbrs(i) ) {// bond_this_atom
			if ( atom_index > i ) {  //don't write out bonds more than once
				out << "BOND  " << rsd.atom_name( i ) << "    " << rsd.atom_name( atom_index ) << " \n";
			}
		}
	} // bond write out


	// now the chis
	for ( Size i=1; i <= rsd.nchi(); ++i ) {
		out << "CHI " << i ;
		AtomIndices atoms_this_chi = rsd.chi_atoms( i );
		for ( AtomIndices::iterator at_it = atoms_this_chi.begin(); at_it != atoms_this_chi.end(); ++at_it ) {
			out << "   " << rsd.atom_name( *at_it );
		}
		out << " \n";
	} //chi write out

	// and now the proton chis
	Size n_proton_chi(0);
	for ( Size i=1; i <= rsd.nchi(); i++ ) {
		if ( rsd.is_proton_chi( i ) ) {

			n_proton_chi++;
			out << "PROTON_CHI " << i << " SAMPLES ";
			utility::vector1< Real > pchi_samples = rsd.proton_chi_samples( n_proton_chi );
			utility::vector1< Real > pchi_extra = rsd.proton_chi_extra_samples( n_proton_chi );
			out << pchi_samples.size() ;
			for ( Size j = 1; j <= pchi_samples.size(); j++ ) { out << " " << pchi_samples[j]; }
			out << " EXTRA " << pchi_extra.size();
			for ( Size j = 1; j <= pchi_extra.size(); j++ ) { out << " " << pchi_extra[j]; }
			out << " \n";

		}
	}//proton chi write out

	// Now the nus...
	for ( Size i = 1, n_nus = rsd.n_nus(); i <= n_nus; ++i ) {
		out << "NU " << i;
		AtomIndices atoms_for_this_nu = rsd.nu_atoms(i);
		for ( AtomIndices::iterator at_it = atoms_for_this_nu.begin(); at_it != atoms_for_this_nu.end(); ++at_it ) {
			out << "   " << rsd.atom_name(*at_it);
		}
		out << std::endl;
	}

	// now all the properties
	out << "PROPERTIES";
	utility::vector1< std::string > const & properties( rsd.properties().get_list_of_properties() );
	Size const n_properties( properties.size() );
	for ( core::uint i = 1; i <= n_properties; ++i ) {
		out << ' ' << properties[ i ];
	}
	out << " \n";

	out << "NBR_ATOM " << rsd.atom_name( rsd.nbr_atom() ) << " \n";
	out << "NBR_RADIUS " << rsd.nbr_radius() << " \n";
	if ( rsd.force_nbr_atom_orient() ) { out << "ORIENT_ATOM NBR\n"; }

	// Charges
	for ( Size i=1; i <= rsd.natoms(); ++i ) {
		if ( rsd.atom(i).formal_charge() != 0 ) {
			out << "CHARGE " << rsd.atom_name( i ) << " FORMAL  " << rsd.atom(i).formal_charge() << " \n";
		}
	}

	// actcoord atoms
	if ( rsd.actcoord_atoms().size() > 0 ) {
		out << "ACT_COORD_ATOMS ";
		AtomIndices act_atoms = rsd.actcoord_atoms();
		for ( AtomIndices::iterator at_it = act_atoms.begin(); at_it != act_atoms.end(); ++at_it ) {
			out << rsd.atom_name( *at_it ) << " ";
		}
		out << "END \n";
	}


	// last but not least the internal coordinates
	for ( Size i=1; i <= rsd.natoms(); i++ ) {
		AtomICoor cur_icoor = rsd.icoor( i );
		out << "ICOOR_INTERNAL   " << rsd.atom_name( i ) << "  " << degrees( cur_icoor.phi() ) << "  ";
		out << degrees( cur_icoor.theta() ) << "  " << cur_icoor.d();
		if ( ( cur_icoor.stub_atom1().atomno() <= rsd.natoms() ) && ( cur_icoor.stub_atom1().atomno() > 0 ) ) {
			out << "   " << rsd.atom_name( cur_icoor.stub_atom1().atomno() );
		} else { out << "   " << cur_icoor.stub_atom1().type(); }

		if ( ( cur_icoor.stub_atom2().atomno() <= rsd.natoms()  ) && ( cur_icoor.stub_atom2().atomno() > 0 ) ) {
			out << "   " << rsd.atom_name( cur_icoor.stub_atom2().atomno() );
		} else { out << "  "  << cur_icoor.stub_atom2().type(); }

		if ( ( cur_icoor.stub_atom3().atomno() <= rsd.natoms()  ) && ( cur_icoor.stub_atom3().atomno() > 0 ) ) {
			out << "   " << rsd.atom_name( cur_icoor.stub_atom3().atomno() );
		} else { out << "  "  << cur_icoor.stub_atom3().type() ; }

		out << " \n";

	} //atom icoor write out

	// TODO: now write out icoors for connections (polymer, other)

	out.close();

} // write_topology_file

/// @brief Callback class for write_graphviz - outputs properties for the nodes and edges.
class GraphvizPropertyWriter {
public:
	GraphvizPropertyWriter( ResidueType const & rsd ):
		rsd_(rsd)
	{}

	/// @brief write properties for the graph as a whole
	void operator()(std::ostream & out) const {
		out << "label=\"" << rsd_.name() << "\"\n";
		//out << "overlap=scale" << "\n";
	}

	/// @brief write properties for atoms
	void operator()(std::ostream & out, VD const & vd) const {
		out << "[";
		Atom const & aa( rsd_.atom(vd) );
		out << "label=\"" << aa.name() << "\"";

		element::Elements elem = aa.element_type()->element();
		if ( elem == element::H ) {
			out << ",color=lightgrey";
		} else if ( elem == element::C ) {
			out << ",color=black";
		} else if ( elem == element::N ) {
			out << ",color=blue";
		} else if ( elem == element::O ) {
			out << ",color=red";
		} else if ( elem == element::S ) {
			out << ",color=gold";
		} else if ( elem == element::P ) {
			out << ",color=orange";
		} else if ( elem == element::F ) {
			out << ",color=chartreuse";
		} else if ( elem == element::Cl ) {
			out << ",color=green";
		} else if ( elem == element::Br ) {
			out << ",color=indianred";
		} else if ( elem == element::I ) {
			out << ",color=indigo";
		} else {
			out << ",color=slategrey";
		}
		out << ",style=bold,shape=circle";
		out << "]";
	}

	/// @brief write properties for bonds
	void operator()(std::ostream & out, ED const & ed) const {
		out << "[";
		Bond const & ee( rsd_.bond(ed) );
		switch ( ee.order() ) {
		case SingleBondOrder : out << "color=black"; break;
		case DoubleBondOrder : out << "color=\"black:white:black\""; break;
		case TripleBondOrder : out << "color=\"black:white:black:white:black\""; break;
		case OrbitalBondOrder : out << "style=dotted"; break;
		case PseudoBondOrder : out << "style=dotted"; break;
		default :
			if ( ee.aromaticity() == IsAromaticBond ) {
				out << "style=dashed";
			} else {
				out << "style=dotted";
			}
		}
		out << "]";
	}

private:
	ResidueType const & rsd_;

};


/// @brief Produces a graphviz dot representation of the ResidueType to the given output stream
/// If header is true (the default) a line with an explanitory message will be printed first.
void
write_graphviz(
	ResidueType const & rsd,
	std::ostream & out,
	bool header /*= true*/
) {
	GraphvizPropertyWriter gpw( rsd );
	if ( header ) {
		out << "// Graphviz dot output for residue " << rsd.name() << std::endl;
		out << "// To use, save to file and run 'neato -T png < residue.dot > residue.png" << "\n"; //supress tracer on next line
	}
	// Go through stringstream to avoid std::endl flushes and tracer name printing.
	std::stringstream ss;
	boost::write_graphviz( ss, rsd.graph(), gpw, gpw, gpw );
	out << ss.str() << std::endl;
} // write_graphviz

////////////////////////////////////////////////////////
void
setup_atom_type_reassignments_from_commandline(
	std::string const & rsd_type_name,
	std::string const & rsd_type_set_name,
	std::map< std::string, std::string > & atom_type_reassignments
)
{
	if ( !basic::options::option[ basic::options::OptionKeys::chemical::reassign_atom_types ].user() ) return;

	if ( rsd_type_name.empty() ) {
		utility_exit_with_message("setup_atom_type_reassignments_from_commandline, empty rsd_type_name");
	}

	utility::vector1< std::string > mods
		( basic::options::option[ basic::options::OptionKeys::chemical::reassign_atom_types ] );

	std::string const errmsg( "-reassign_atom_types format should be:: -reassign_atom_types <rsd-type-set1-name>:<rsd-type1-name>:<atom1-name>:<new-atom-type1-name>   <rsd-type-set2-name>:<rsd-type2-name>:<atom2-name>:<new-atom-type2-name> ...; for example: '-chemical:reassign_atom_types fa_standard:ARG:NE:NtpR' ");

	for ( Size i=1; i<= mods.size(); ++i ) {
		///
		/// mod should look like (for example):  "fa_standard:OOC:LK_RADIUS:4.5"
		///
		std::string const & mod( mods[i] );

		Size const pos1( mod.find(":") );
		if ( pos1 == std::string::npos ) utility_exit_with_message(errmsg);
		std::string const mod_rsd_type_set_name( mod.substr(0,pos1) );
		if ( mod_rsd_type_set_name != rsd_type_set_name ) continue;

		Size const pos2( mod.substr(pos1+1).find(":") );
		if ( pos2 == std::string::npos ) utility_exit_with_message(errmsg);
		std::string const mod_rsd_type_name( mod.substr(pos1+1,pos2) );
		if ( mod_rsd_type_name != rsd_type_name ) continue;

		Size const pos3( mod.substr(pos1+1+pos2+1).find(":") );
		if ( pos3 == std::string::npos ) utility_exit_with_message(errmsg);
		std::string const atom_name( mod.substr(pos1+1+pos2+1,pos3) );

		std::string const new_atom_type_name( mod.substr(pos1+1+pos2+1+pos3+1) );

		tr.Trace << "setup_atom_type_reassignments_from_commandline: reassigning " << rsd_type_set_name << ' ' << rsd_type_name << ' ' << atom_name << " to new atomtype: " << new_atom_type_name << std::endl;

		atom_type_reassignments[ atom_name ] = new_atom_type_name;
	}


}
////////////////////////////////////////////////////////
void
setup_atomic_charge_reassignments_from_commandline(
	std::string const & rsd_type_name,
	std::string const & rsd_type_set_name,
	std::map< std::string, Real > & atomic_charge_reassignments
)
{
	if ( !basic::options::option[ basic::options::OptionKeys::chemical::set_atomic_charge ].user() ) return;

	if ( rsd_type_name.empty() ) {
		utility_exit_with_message("setup_atomic_charge_reassignments_from_commandline, empty rsd_type_name");
	}

	utility::vector1< std::string > mods
		( basic::options::option[ basic::options::OptionKeys::chemical::set_atomic_charge ] );

	std::string const errmsg( "-set_atomic_charge format should be::  -chemical:set_atomic_charge <rsd-type-set1-name>:<rsd-type1-name>:<atom1-name>:<new-charge> <rsd-type-set2-name>:<rsd-type2-name>:<atom2-name>:<new-charge>  ... For example: '-chemical:set_atomic_charge fa_standard:ARG:NE:-1' ");

	for ( Size i=1; i<= mods.size(); ++i ) {
		///
		/// mod should look like (for example):  "fa_standard:OOC:LK_RADIUS:4.5"
		///
		std::string const & mod( mods[i] );

		Size const pos1( mod.find(":") );
		if ( pos1 == std::string::npos ) utility_exit_with_message(errmsg);
		std::string const mod_rsd_type_set_name( mod.substr(0,pos1) );
		if ( mod_rsd_type_set_name != rsd_type_set_name ) continue;

		Size const pos2( mod.substr(pos1+1).find(":") );
		if ( pos2 == std::string::npos ) utility_exit_with_message(errmsg);
		std::string const mod_rsd_type_name( mod.substr(pos1+1,pos2) );
		if ( mod_rsd_type_name != rsd_type_name ) continue;

		Size const pos3( mod.substr(pos1+1+pos2+1).find(":") );
		if ( pos3 == std::string::npos ) utility_exit_with_message(errmsg);
		std::string const atom_name( mod.substr(pos1+1+pos2+1,pos3) );

		std::string const new_atomic_charge_string( mod.substr(pos1+1+pos2+1+pos3+1) );

		if ( !ObjexxFCL::is_float( new_atomic_charge_string ) ) utility_exit_with_message(errmsg);
		Real const new_atomic_charge( ObjexxFCL::float_of( new_atomic_charge_string ) );

		tr.Trace << "setup_atomic_charge_reassignments_from_commandline: setting charge of " << rsd_type_set_name << ' ' << rsd_type_name << ' ' << atom_name << " to " << new_atomic_charge << std::endl;

		atomic_charge_reassignments[ atom_name ] = new_atomic_charge;
	}


}
////////////////////////////////////////////////////////
void
setup_icoor_reassignments_from_commandline(
	std::string const & rsd_type_name,
	std::string const & rsd_type_set_name,
	std::map< std::string, utility::vector1< std::string > > & icoor_reassignments
)
{
	if ( !basic::options::option[ basic::options::OptionKeys::chemical::reassign_icoor ].user() ) return;

	if ( rsd_type_name.empty() ) {
		utility_exit_with_message("setup_icoor_reassignments_from_commandline, empty rsd_type_name");
	}

	utility::vector1< std::string > mods
		( basic::options::option[ basic::options::OptionKeys::chemical::reassign_icoor ] );

	std::string const errmsg( "-reassign_icoor format should be:: -reassign_icoor <rsd-type-set1-name>:<rsd-type1-name>:<atom1-name>:<the-six-icoor-params-as-a-comma-separated-list>   <rsd-type-set2-name>:<rsd-type2-name>:<atom2-name>:<icoorparams2> ...; for example: -chemical:reassign_icoor fa_standard:ADE:UPPER:-180,60,1.6,O3',C3',C4' ");

	for ( Size i=1; i<= mods.size(); ++i ) {
		///
		///
		std::string const mod( stripped( mods[i] ) );

		Size const pos1( mod.find(":") );
		if ( pos1 == std::string::npos ) utility_exit_with_message(errmsg);
		std::string const mod_rsd_type_set_name( mod.substr(0,pos1) );
		if ( mod_rsd_type_set_name != rsd_type_set_name ) continue;

		Size const pos2( mod.substr(pos1+1).find(":") );
		if ( pos2 == std::string::npos ) utility_exit_with_message(errmsg);
		std::string const mod_rsd_type_name( mod.substr(pos1+1,pos2) );
		if ( mod_rsd_type_name != rsd_type_name ) continue;

		Size const pos3( mod.substr(pos1+1+pos2+1).find(":") );
		if ( pos3 == std::string::npos ) utility_exit_with_message(errmsg);
		std::string const atom_name( mod.substr(pos1+1+pos2+1,pos3) );

		utility::vector1< std::string > const new_icoor_params(utility::string_split( mod.substr(pos1+1+pos2+1+pos3+1),','));
		runtime_assert( new_icoor_params.size() == 6 );
		runtime_assert( ObjexxFCL::is_float( new_icoor_params[1] ) );
		runtime_assert( ObjexxFCL::is_float( new_icoor_params[2] ) );
		runtime_assert( ObjexxFCL::is_float( new_icoor_params[3] ) );

		tr.Trace << "setup_icoor_reassignments_from_commandline: reassigning " << rsd_type_set_name << ' ' <<
			rsd_type_name << ' ' << atom_name << " to new icoor params: ";
		for ( Size i=1; i<= 6; ++i ) tr.Trace << ' ' << new_icoor_params[i];
		tr.Trace << std::endl;

		icoor_reassignments[ atom_name ] = new_icoor_params;
	}
}

/// @brief Symmetrize the glycine params file (if the user has used the -symmetric_gly_tables option).
/// @details Ugh.  Special-case logic.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
void
apply_symm_gly_corrections(
	std::string const &child_atom,
	core::Real &phi,
	core::Real &/*theta*/,
	core::Real &d,
	std::string &/*parent_atom*/,
	std::string &/*angle_atom*/,
	std::string &torsion_atom
) {

	if ( child_atom == "UPPER" ) {
		phi = 180.0;
	} else if ( child_atom == "1HA" ) {
		d = 1.089761;
	} else if ( child_atom == "2HA" ) {
		phi = -121.4;
		d = 1.089761;
		torsion_atom = "C";
	} else if ( child_atom == "LOWER" ) {
		phi = 180.0;
	}

	return;
}

} // chemical
} // core
