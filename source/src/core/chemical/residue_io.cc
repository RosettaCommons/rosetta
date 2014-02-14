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
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSupport.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/MMAtomType.hh>

// Project headers
#include <platform/types.hh>

#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/Bound.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/SmallKeyVector.hh>

// External headers
#include <boost/foreach.hpp>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray2D.hh>


#define foreach BOOST_FOREACH

namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS


namespace core {
namespace chemical {

static basic::Tracer tr("core.chemical");


///////////////////////////////////////////////////////////////////////////////
/// @brief helper fxn
id::AtomID
atom_id_from_icoor_line(
	std::string const name,
	ResidueType const & rsd
)
{
	using id::AtomID;
	ICoorAtomID id( name, rsd );

	switch ( id.type() ) {
	case ICoorAtomID::INTERNAL:
		return AtomID( id.atomno(), 1 );
	case ICoorAtomID::CONNECT:
		return AtomID( id.atomno(), 2 );
	case ICoorAtomID::POLYMER_LOWER:
		return AtomID( 1, 3 );
	case ICoorAtomID::POLYMER_UPPER:
		return AtomID( 2, 3 );
	default:
		utility_exit_with_message( "unrecognized stub atom id type!" );
	}
	return id::BOGUS_ATOM_ID;
}


ResidueTypeOP
read_topology_file(
	std::string const & filename,
	chemical::AtomTypeSetCAP atom_types,
	chemical::ElementSetCAP elements,
	chemical::MMAtomTypeSetCAP mm_atom_types,
	chemical::orbitals::OrbitalTypeSetCAP orbital_atom_types,
//	chemical::CSDAtomTypeSetCAP csd_atom_types kwk commenting out csd_atom_types until I have a chance to fully implement them.
	chemical::ResidueTypeSetCAP rsd_type_set
)
{
	if( ! utility::file::file_exists( filename ) ) {
		utility_exit_with_message("Cannot find file '"+filename+"'");
	}
	utility::io::izstream data( filename.c_str() );
	if ( !data.good() ) {
		utility_exit_with_message("Cannot open file '"+filename+"'");
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
/// BOND:
/// Declares a bond between two atoms giving their names. This line is
/// whitespace delimited.  E.g., "BOND  N    CA" from ALA.params.
///
/// BOND_TYPE:
/// Declares a bond between two atoms, giving their names, and also
/// describing the chemical nature of the bond.  See the "BondName"
/// enumeration for acceptable bond types.
///
/// CHARGE:
/// Declares a charge for a particular atom.
/// Format CHARGE atom type value
/// Currently valid types are FORMAL. (Partial charges are handled above.)
/// E.g. "CHARGE OG2 FORMAL -1"
///
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
/// should be listed before any additional CONNECT records.
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
/// are considered backbone atoms.  All heavy atoms after the first
/// side-chain atom are considered side-chain atoms. Hydrogen atoms
/// are either side-chain or backbone depending on the heavy atom to
/// which they are bound. E.g., "FIRST_SIDECHAIN_ATOM CB" from SER.params.
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
/// connections, atoms "CONN#" should be given (e.g. CONN3 for the disulfide connection in CYD).
/// The number for an inter-residue connection comes from the order in which the connection is
/// declared in the file, and includes the LOWER_CONNECT and UPPER_CONNECT connections in this
/// count (e.g. for CYD, there is a LOWER_CONNECT, and UPPER_CONNECT, and only a single
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
/// NAME:
/// Gives the name for this ResidueType.  The name for each ResidueType
/// must be unique within a ResidueTypeSet.  It is not limited to three latters.
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
///
/// NCAA_ROTLIB_NUM_ROTAMER_BINS:
/// Lists the number of rotamers and the number of bins for each rotamer.
/// E.g., "NCAA_ROTLIB_NUM_ROTAMER_BINS 2 3 2" from
/// d-ncaa/d-5-fluoro-tryptophan.params
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
/// a residue.  This is commonly used for small molecules.  See also
/// the CHI_ROTAMERS and NCAA_ROTLIB_PATH tags for alternate ways to
/// provide rotamers.
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
/// ROTAMER_AA:
/// Sets the "rotamer_aa" for a particular residue, which can be used
/// to describe to the RotamerLibrary what amino acid to mimic for the
/// sake of building rotamers.  E.g., "ROTAMER_AA SER" No examples
/// currently found in the database (10/13).
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
/// other. E.g. "VIRTUAL_SHADOW NV N" from PRO.params.
ResidueTypeOP
read_topology_file(
		utility::io::izstream & data,
		chemical::AtomTypeSetCAP atom_types,
		chemical::ElementSetCAP elements,
		chemical::MMAtomTypeSetCAP mm_atom_types,
		chemical::orbitals::OrbitalTypeSetCAP orbital_atom_types,
		//chemical::CSDAtomTypeSetCAP csd_atom_types kwk commenting out until they have been fully implemented
		chemical::ResidueTypeSetCAP rsd_type_set
)
{
	assert( rsd_type_set );

	using id::AtomID;
	using id::DOF_ID;
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	using namespace basic;

	std::string filename = data.filename();

	// read the file
    std::string line;
	utility::vector1< std::string > lines;

	while ( getline( data, line ) ) {
		std::istringstream l( line );
		//if ( line.size() < 1 || line[0] == '#' ) continue;
		if ( line.size() < 1 ) continue;
		std::string::size_type pound = line.find('#', 0);
		if( pound == std::string::npos ) {
			lines.push_back( line );
		} else {
			std::string no_comment_line= line.substr(0, pound);
			lines.push_back(no_comment_line);
		}
	}
	tr.Debug << "Read " << lines.size() << " lines from file: " << filename << std::endl;
	data.close();


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

	ResidueTypeOP rsd( new ResidueType( atom_types, elements, mm_atom_types, orbital_atom_types ) ); //kwk commenting out until atom types are fully implemented , csd_atom_types ) );
	rsd->residue_type_set( rsd_type_set );  // Give this rsd_type a backpointer to its set.

	// Add the atoms.
	Size const nlines( lines.size() );
	Size natoms(0);//, norbitals(0);
	for (Size i=1; i<= nlines; ++i) {
		std::string line( lines[i] );
		if (line.size() > 0) {
			while (line.substr(line.size()-1) == " ") {
				line = line.substr(0, line.size()-1);
				if (line.size() == 0) break;
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
		l >> tag; std::string const atom_type_name( tag );

		// read in the Molecular mechanics atom type
		// std::string const mm_atom_type_name( line.substr(15,4) );
		l >> tag; std::string const mm_atom_type_name( tag );

		// the atomic charge
		float charge;
		// std::istringstream l( line.substr(20) );
		l >> charge;
		float parse_charge(charge);
		if (!l.eof()) {
			l >> parse_charge;
		}

		if ( ! basic::options::option[ basic::options::OptionKeys::corrections::chemical::parse_charge ]() ) {
			rsd->add_atom( atom_name, atom_type_name, mm_atom_type_name, charge );
		}
		else {
			rsd->add_atom( atom_name, atom_type_name, mm_atom_type_name, parse_charge );
		}

		++natoms;
	}

	// No ATOM lines probably means an invalid file.  Perhaps someone made a mistake with an -extra_res_fa flag.
	// Fail gracefully now, versus a segfault later.
	if (natoms == 0) {
		utility_exit_with_message("Residue topology file '" + filename + "' does not contain valid ATOM records.");
	}

	// Add the bonds; parse the rest of file.
	bool found_AA_record = false;
	bool found_PDB_ROTAMERS_record = false;
	std::string pdb_rotamers_filename = "";

	for (Size i=1; i<= nlines; ++i) {
		std::string const & line( lines[i] );
		std::istringstream l( line );
		std::string tag,atom1,atom2,atom3,atom4, rotate, orbitals_tag, orbital;
		core::Real value;
		core::Size bond_type;
		l >> tag;
		if ( l.fail() ) continue;
		if ( tag == "CONNECT" ) {
			l >> atom1;
			l >> rotate; // not used here
			rsd->add_residue_connection( atom1);
			//std::cout << "CONNECT record deprecated " << std::endl;
		} else if ( tag == "TYPE" ) {
			// will probably handle this differently later on
			l >> tag;
			if ( tag == "POLYMER" ) {
				rsd->add_property( tag );
			} else if ( tag == "LIGAND" ) {
				rsd->add_property( tag );
			}
		} else if ( tag == "BOND" ) {
			l >> atom1 >> atom2;
			rsd->add_bond( atom1, atom2 );

		} else if ( tag == "BOND_TYPE" ) {
			l >> atom1 >> atom2 >> bond_type;
			// apl Note: this cast could easily fail and there's no error checking
			rsd->add_bond(atom1, atom2, static_cast<core::chemical::BondName>(bond_type));
		} else if ( tag == "CHARGE" ) {
			l >> atom1;
			// We should allow multiple charges on one line, but since we now just have the one, hold off.
			l >> tag;
			if( tag == "FORMAL" ) {
				l >> value;
				rsd->atom( atom1 ).formal_charge( int(value) );
			} else {
				utility_exit_with_message("ERROR: Invalid charge type '"+tag+"' in topology file.");
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
		} else if ( tag == "PROTON_CHI") {
			Size chino, nsamples, nextra_samples;
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
			if ( basic::options::option[ basic::options::OptionKeys::corrections::chemical::expand_st_chi2sampling ]
					&& (rsd->aa() == aa_ser || rsd->aa() == aa_thr )
					&& rsd_type_set->name() == FA_STANDARD ) {
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
				utility_exit_with_message("Unknown ORIENT_ATOM mode: " + tag );
			}
		} else if ( tag == "PROPERTIES" ) {
			l >> tag;
			while ( !l.fail() ) {
				rsd->add_property( tag );
				l >> tag;
			}
		} else if (tag == "NUMERIC_PROPERTY"){
			core::Real value = 0.0;
			l >> tag >> value;
			rsd->add_numeric_property(tag,value);

		} else if (tag == "STRING_PROPERTY" ) {
			std::string value;
			l >> tag >> value;
			rsd->add_string_property(tag,value);

		} else if ( tag == "VARIANT" ) {
			l >> tag;
			while ( !l.fail() ) {
				rsd->add_variant_type( tag );
				l >> tag;
			}
		} else if ( tag == "FIRST_SIDECHAIN_ATOM" ) {
			// note-- atoms are sidechain by default
			l >> tag;
			if ( tag == "NONE" ) {
				// set all atoms to backbone
				for ( Size j=1; j<= rsd->natoms(); ++j ) {
					rsd->set_backbone_heavyatom( rsd->atom_name(j) );
				}
			} else if ( rsd->has( tag ) ) {
				for ( Size j=1; j< rsd->atom_index( tag ); ++j ) {
					rsd->set_backbone_heavyatom( rsd->atom_name(j) );
				}
			}

		} else if ( tag == "IO_STRING" ) {
			assert( line.size() >= 15 );
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

		} else if ( tag == "ROTAMER_AA" ) {
			l >> tag;
			rsd->rotamer_aa( tag );

		} else if ( tag == "NAME" ) {
			l >> tag;
			rsd->name( tag );

		} else if ( tag == "CHI_ROTAMERS" ) {
			Size chino;
			Real mean, sdev;
			l >> chino;
			l >> mean >> sdev;
			while ( !l.fail() ) {
				rsd->add_chi_rotamer( chino, mean, sdev );
				l >> mean >> sdev;
			}
		} else if ( tag == "PDB_ROTAMERS" ) {
			found_PDB_ROTAMERS_record = true;
			l >> pdb_rotamers_filename;
		} else if ( tag == "ACT_COORD_ATOMS" ) {
			while ( l ) {
				l >> atom1;
				if ( atom1 == "END") break;
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
		} else if ( tag == "NCAA_ROTLIB_PATH" ) {
			std::string path;
			l >> path;
			rsd->set_ncaa_rotlib_path( path );
			rsd->set_use_ncaa_rotlib( true );
		} else if ( tag == "NCAA_ROTLIB_NUM_ROTAMER_BINS" ) {
			Size n_rots(0);
			utility::vector1<Size> n_bins_per_rot;
			l >> n_rots;
			rsd->set_ncaa_rotlib_n_rotameric_bins( n_rots );
			n_bins_per_rot.resize( n_rots );
			for( Size i = 1; i <= n_rots; ++i ) {
				Size bin_size(0);
				l >> bin_size;
				n_bins_per_rot[i] = bin_size;
			}
			rsd->set_ncaa_rotlib_n_bin_per_rot( n_bins_per_rot );
		}  else if ( tag == "VIRTUAL_SHADOW" ) {
			std::string shadower, shadowee;
			l >> shadower >> shadowee;
			rsd->set_shadowing_atom( shadower, shadowee );
		} else if ( tag == "ATOM" || tag == "ICOOR_INTERNAL" ){
			; // ATOM lines handled above, ICOOR_INTERNAL lines handled below
		} else {
			tr.Warning << "WARNING: Ignoring line starting with '" << tag << "' when parsing topology file." << std::endl;
		}

	} // i=1,nlines


	if ( !found_AA_record ) {
		basic::Warning() << "No AA record found for " << rsd->name()
				<< "; assuming " << name_from_aa( rsd->aa() ) << std::endl;
	}


	// set icoor coordinates, store information about polymer links
	// also sets up base_atom
	{

		std::map< std::string, Vector > rsd_xyz;

		for ( Size i=1; i<= nlines; ++i ) {

			std::string const & line( lines[i] );
			std::istringstream l( line );
			std::string tag, child_atom, parent_atom, angle_atom, torsion_atom;

			Real phi, theta, d;
			l >> tag;

			if ( tag != "ICOOR_INTERNAL" ) continue;

			l >> child_atom >> phi >> theta >> d >> parent_atom >> angle_atom >> torsion_atom;

			phi = radians(phi); theta = radians(theta); // in degrees in the file for human readability

			// This code should probably be extracted to a util function
			if ( natoms > 1 ) {
				// build the Cartesian coords for the new atom:
				if ( child_atom == parent_atom ) {
					if( ! rsd_xyz.empty() ) {
						utility_exit_with_message("Only the first ICOOR atom in a topology file should list itself as its own parent atom.");
					}
					rsd_xyz[ child_atom ] = Vector( 0.0 );

				} else if ( child_atom == angle_atom ) {
					if( rsd_xyz.size() != 1 ) {
						utility_exit_with_message("Only the second ICOOR atom in a topology file should list itself as its own angle atom.");
					}
					if( ! rsd_xyz.count( parent_atom ) ) {
						utility_exit_with_message("In second ICOOR atom in topology file - parent atom not found.");
					}
					rsd_xyz[ child_atom ] = Vector( d, 0.0, 0.0 );

				} else {
					Vector torsion_xyz;
					if ( child_atom == torsion_atom ) {
						if( rsd_xyz.size() != 2 ) {
							utility_exit_with_message("Only the third ICOOR atom in a topology file should list itself as its own dihedral atom.");
						}
						if( ! rsd_xyz.count( parent_atom ) || ! rsd_xyz.count( angle_atom ) ) {
							utility_exit_with_message("In third ICOOR atom in topology file - parent and/or angle atom not found.");
						}
						torsion_xyz = Vector( 1.0, 1.0, 0.0 );
					} else {
						if( ! ( rsd_xyz.count( parent_atom ) && rsd_xyz.count( angle_atom ) &&
								rsd_xyz.count( torsion_atom ) ) ) {
							utility_exit_with_message("In ICOOR atom line in topology file: reference atoms must be specified in earlier line. One of "+parent_atom+" or "+angle_atom+" or "+torsion_atom);
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
			rsd->set_icoor( i, child_atom, phi, theta, d, parent_atom, angle_atom, torsion_atom );


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
				assert( rsd_xyz.count( name ) );
				rsd->set_ideal_xyz( name, rsd_xyz[ name ] );
				//rsd->set_xyz( rsd->atom_name(i), atom_tree.xyz( id::AtomID(i,1) ) );
				//rsd->atom(i).xyz( atom_tree.xyz( id::AtomID(i,1) ) );
			}
		}

		// If polymer, fill list of main chain atoms.
		if ( rsd->is_polymer() ) {
			// Test that this is really a polymer residue.
			if ( rsd->upper_connect_id() && rsd->upper_connect_atom() &&
					rsd->lower_connect_id() && rsd->lower_connect_atom() ) {
				Size upper_connect( rsd->upper_connect_atom() ), lower_connect( rsd->lower_connect_atom() );
				AtomIndices mainchain;
				if (!rsd->is_NA()) {  // Default main chain: defined by order in .params file
					for (uint atom = lower_connect; atom <= upper_connect; ++atom) {
						mainchain.push_back(atom);
					}
				} else /* rsd->is_NA() */ {  // Alternative main chain: defined by shortest path from LOWER to UPPER
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
						assert( next_atom != atom );

						atom = next_atom;
					}
					mainchain.push_back( upper_connect );
				}
				rsd->set_mainchain_atoms( mainchain );
			} else {
				tr.Warning << "WARNING: Residue " << rsd->name() << " claims it's a polymer, " <<
						"but it doesn't have the appropriate UPPER and LOWER connection points specified." << std::endl;
			}
		}

		// now also need to store the information about the geometry at the links...

	} // scope

	// calculate any remaining derived data
	rsd->finalize();

	// If we found a PDB_ROTAMERS library, load them now that the ResidueType
	// is totally initialized:
	if( found_PDB_ROTAMERS_record ) {
		using namespace utility::file;
		// Assume name of rotamers file has no path info, etc
		// and is in same directory as the residue parameters file.
		FileName this_file( filename ), rot_file( pdb_rotamers_filename );
		rot_file.vol( this_file.vol() );
		rot_file.path( this_file.path() );

		rsd->set_RotamerLibraryName( rot_file() );

		tr.Debug << "Setting up conformer library for " << rsd->name() << std::endl;
		/*using namespace core::pack::dunbrack;
		using namespace utility::file;
		SingleLigandRotamerLibraryOP pdb_rotamers = new SingleLigandRotamerLibrary();
		// Assume name of rotamers file has no path info, etc
		// and is in same directory as the residue parameters file.

		pdb_rotamers->init_from_file( rot_file.name(), rsd );
		rsd->set_RotamerLibrary( pdb_rotamers );*/
	}

	return rsd;
}


/// @brief function to write out a topology file given a residue type, can be used to
/// @brief debug on the fly generated residue types. Note: not perfect yet, the enums for
/// @brief the connection types are given in numbers instead of names
void
write_topology_file(
	ResidueType const & rsd
)
{

	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	std::string filename = rsd.name() + ".params";

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
	for (Size i=1; i <= rsd.natoms(); ++i){

		std::string atom_out = "ATOM " + rsd.atom_name( i ) + " " + rsd.atom_type( i ).name() + "  ";
		atom_out = atom_out + rsd.mm_atom_type(i).name();
		out << atom_out << " " << rsd.atom(i).charge() << " \n";

	} // atom write out

	if ( rsd.is_polymer() ) {
		if ( !rsd.is_lower_terminus() ) { out << "LOWER_CONNECT " << rsd.atom_name( rsd.lower_connect().atomno() ) << " \n"; }
		if ( !rsd.is_upper_terminus() ) { out << "UPPER_CONNECT " << rsd.atom_name( rsd.upper_connect().atomno() ) << " \n";}
	}

	// then all the bonds
	for (Size i=1; i <= rsd.natoms(); ++i){
		foreach(Size atom_index, rsd.nbrs(i)){// bond_this_atom
			if( atom_index > i ) {  //don't write out bonds more than once
				out << "BOND  " << rsd.atom_name( i ) << "    " << rsd.atom_name( atom_index ) << " \n";
			}
		}
	} // bond write out


	// now the chis
	for (Size i=1; i <= rsd.nchi(); ++i){
		out << "CHI " << i ;
		AtomIndices atoms_this_chi = rsd.chi_atoms( i );
		for (AtomIndices::iterator at_it = atoms_this_chi.begin(); at_it != atoms_this_chi.end(); ++at_it){
			out << "   " << rsd.atom_name( *at_it );
		}
		out << " \n";
	} //chi write out

	// and now the proton chis
	Size n_proton_chi(0);
	for(Size i=1; i <= rsd.nchi(); i++){
		if( rsd.is_proton_chi( i ) ){

			n_proton_chi++;
			out << "PROTON_CHI " << i << " SAMPLES ";
			utility::vector1< Real > pchi_samples = rsd.proton_chi_samples( n_proton_chi );
			utility::vector1< Real > pchi_extra = rsd.proton_chi_extra_samples( n_proton_chi );
			out << pchi_samples.size() ;
			for( Size j = 1; j <= pchi_samples.size(); j++) { out << " " << pchi_samples[j]; }
			out << " EXTRA " << pchi_extra.size();
			for( Size j = 1; j <= pchi_extra.size(); j++) { out << " " << pchi_extra[j]; }
			out << " \n";

		}
	}//proton chi write out

	// Now the nus...
	for (Size i = 1, n_nus = rsd.n_nus(); i <= n_nus; ++i){
		out << "NU " << i;
		AtomIndices atoms_for_this_nu = rsd.nu_atoms(i);
		for (AtomIndices::iterator at_it = atoms_for_this_nu.begin(); at_it != atoms_for_this_nu.end(); ++at_it) {
			out << "   " << rsd.atom_name(*at_it);
		}
		out << std::endl;
	}

	// now all the properties
	out << "PROPERTIES";
	if (rsd.is_protein() ) { out << " PROTEIN"; }
	if (rsd.is_DNA() ) { out << " DNA"; }
	if (rsd.is_RNA() ) { out << " RNA"; }
	if (rsd.is_polar() ) { out << " POLAR"; }
	if (rsd.is_charged() ) { out << " CHARGED"; }
	if (rsd.is_aromatic() ) { out << " AROMATIC"; }
	if (rsd.is_lower_terminus() ) { out << " LOWER_TERMINUS"; }
	if (rsd.is_upper_terminus() ) { out << " UPPER_TERMINUS"; }
	if (rsd.is_terminus() ) { out << " TERMINUS"; }
	out << " \n";

	out << "NBR_ATOM " << rsd.atom_name( rsd.nbr_atom() ) << " \n";
	out << "NBR_RADIUS " << rsd.nbr_radius() << " \n";
	if (rsd.force_nbr_atom_orient()) { out << "ORIENT_ATOM NBR\n"; }

	// Charges
	for (Size i=1; i <= rsd.natoms(); ++i){
		if( rsd.atom(i).formal_charge() != 0 ) {
			out << "CHARGE " << rsd.atom_name( i ) << " FORMAL  " << rsd.atom(i).formal_charge() << " \n";
		}
	}

	// actcoord atoms
	if ( rsd.actcoord_atoms().size() > 0 ){
		out << "ACT_COORD_ATOMS ";
		AtomIndices act_atoms = rsd.actcoord_atoms();
		for(AtomIndices::iterator at_it = act_atoms.begin(); at_it != act_atoms.end(); at_it++ ){
			out << rsd.atom_name( *at_it ) << " ";
		}
		out << "END \n";
	}


	// last but not least the internal coordinates
	for (Size i=1; i <= rsd.natoms(); i++){
		AtomICoor cur_icoor = rsd.icoor( i );
		out << "ICOOR_INTERNAL   " << rsd.atom_name( i ) << "  " << degrees( cur_icoor.phi() ) << "  ";
		out << degrees( cur_icoor.theta() ) << "  " << cur_icoor.d();
		if( ( cur_icoor.stub_atom1().atomno() <= rsd.natoms() ) && ( cur_icoor.stub_atom1().atomno() > 0 ) ) {
			out << "   " << rsd.atom_name( cur_icoor.stub_atom1().atomno() );
		}
		else{ out << "   " << cur_icoor.stub_atom1().type(); }

		if( ( cur_icoor.stub_atom2().atomno() <= rsd.natoms()  ) && ( cur_icoor.stub_atom2().atomno() > 0 )){
			out << "   " << rsd.atom_name( cur_icoor.stub_atom2().atomno() );
		}
		else{ out << "  "  << cur_icoor.stub_atom2().type();}

		if( ( cur_icoor.stub_atom3().atomno() <= rsd.natoms()  ) && ( cur_icoor.stub_atom3().atomno() > 0 )){
			out << "   " << rsd.atom_name( cur_icoor.stub_atom3().atomno() );
		}
		else{ out << "  "  << cur_icoor.stub_atom3().type() ;}

		out << " \n";

	} //atom icoor write out

	// TODO: now write out icoors for connections (polymer, other)

	out.close();

} // write_topology_file

} // chemical
} // core
