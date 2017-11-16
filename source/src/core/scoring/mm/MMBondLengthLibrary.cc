// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/mm/MMBondLengthLibary.cc
/// @brief  Molecular mechanics bond length score class
/// @author Frank DiMaio (based on Colin Smith's MMBondAngle potential)

// Unit headers
#include <core/scoring/mm/MMBondLengthLibrary.hh>

// Project headers
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers

// C++ headers
#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include <utility/assert.hh>


#include <core/chemical/MMAtomType.hh>


namespace core {
namespace scoring {
namespace mm {

/// @details Auto-generated virtual destructor
MMBondLengthLibrary::~MMBondLengthLibrary() {}

using basic::Error;
using basic::Warning;
static basic::Tracer TR( "core.mm.MMBondLengthLibrary" );

/// @details Construct a MMBondLengthLibrary instant from a filename string and constant access pointer to an MMAtomTypeSet
MMBondLengthLibrary::MMBondLengthLibrary(
	std::string filename,
	core::chemical::MMAtomTypeSetCOP mm_atom_set
)
{
	mm_atom_set_ = core::chemical::MMAtomTypeSetCAP( mm_atom_set );

	// read the file
	std::string line;
	utility::vector1< std::string > lines;
	std::ifstream data( filename.c_str() );

	bool in_bonds_section = false;
	while ( getline( data, line ) ) {
		if ( line.size() < 1 || line[0] == '!' || line[0] == ' ' ) continue; // comment or blank lines
		if ( line == "ANGLES" ) in_bonds_section = false;
		if ( in_bonds_section ) lines.push_back( line );
		if ( line == "BONDS" ) in_bonds_section = true;
	}

	// parse params
	for ( Size i = 1; i <= lines.size(); ++i ) {
		std::istringstream l( lines[i] );

		std::string atom_type_string_1, atom_type_string_2;
		l >> atom_type_string_1 >> atom_type_string_2;

		// skip the parameters if any of the atom types don't exist
		if ( ! mm_atom_set->contains_atom_type( atom_type_string_1 ) ) continue;
		if ( ! mm_atom_set->contains_atom_type( atom_type_string_2 ) ) continue;

		// get atom-type_index from atom set
		int atom_type_int1 = mm_atom_set->atom_type_index( atom_type_string_1 );
		int atom_type_int2 = mm_atom_set->atom_type_index( atom_type_string_2 );

		// get k_b and b_0
		Real k_b, b_0;
		l >> k_b >> b_0;

		mm_bondlength_library_.insert( std::make_pair(
			mm_bondlength_atom_pair( atom_type_int1, atom_type_int2 ),
			mm_bondlength_param_set( k_b, b_0 ) ) );
	}

	/// apl -- add "no-op" pair for virtual atoms
	int const virt_type = mm_atom_set->atom_type_index("VIRT");
	Real const noop_kb( 0.0 );
	Real const noop_b0( 0.0 );
	mm_bondlength_library_.insert( std::make_pair(
		mm_bondlength_atom_pair( virt_type, virt_type ),
		mm_bondlength_param_set( noop_kb, noop_b0 ) ));


	// print number torsion params added
	TR << "MM bond length sets added: " << mm_bondlength_library_.size()
		<< " (+1 virtual)." << std::endl;
}

mm_bondlength_library_citer_pair
MMBondLengthLibrary::lookup( int atom1, int atom2 ) const {
	static std::string const virt_string = "VIRT";

	if ( mm_bondlength_library_.count( mm_bondlength_atom_pair( atom1, atom2 ) ) ) {
		// forward
		return mm_bondlength_library_.equal_range( mm_bondlength_atom_pair( atom1, atom2 ) );
	} else if ( mm_bondlength_library_.count( mm_bondlength_atom_pair( atom2, atom1 ) ) ) {
		// backward
		return mm_bondlength_library_.equal_range( mm_bondlength_atom_pair( atom2, atom1 ) );
	}
	core::chemical::MMAtomTypeSetCOP mm_atom_set( mm_atom_set_ );
	int const virt_atom_type = mm_atom_set->atom_type_index( virt_string );

	/// Virtual atoms get no mm-parameters.  Return the no-op torsion object
	if ( atom1 == virt_atom_type || atom2 == virt_atom_type ) {
		return mm_bondlength_library_.equal_range( mm_bondlength_atom_pair( virt_atom_type, virt_atom_type ));
	}

	TR << "No parameters for " << (*mm_atom_set)[atom1].name() << "-" << (*mm_atom_set)[atom2].name() << std::endl;
	utility_exit_with_message("COULD NOT FIND BOND LENGTH PARAMS" );

	return mm_bondlength_library_citer_pair(); ///< meaningless, just for removing gcc warning.
}

mm_bondlength_library_citer_pair
MMBondLengthLibrary::lookup ( std::string atom1, std::string atom2 ) const {
	core::chemical::MMAtomTypeSetCOP mm_atom_set( mm_atom_set_ );
	return (*this).lookup(
		mm_atom_set->atom_type_index( atom1 ),
		mm_atom_set->atom_type_index( atom2 ) );
}

void
MMBondLengthLibrary::pretty_print() const {
	// for each key print out its value
	for ( mm_bondlength_library_citer i = mm_bondlength_library_.begin(),
			e = mm_bondlength_library_.end(); i != e; ++i ) {
		TR << (i->first).key1() << "\t"
			<< (i->first).key2() << "\t"
			<< (i->second).key1() << "\t"
			<< (i->second).key2() << "\t"
			<< std::endl;
	}
}

void
MMBondLengthLibrary::pretty_print( int atom1, int atom2 ) const {
	mm_bondlength_library_citer_pair temppair = this->lookup(atom1, atom2);
	for ( mm_bondlength_library_citer i = temppair.first, e = temppair.second; i != e; ++i ) {
		TR << (i->first).key1() << "\t"
			<< (i->first).key2() << "\t"
			<< (i->second).key1() << "\t"
			<< (i->second).key2() << "\t"
			<< std::endl;
	}
}

void
MMBondLengthLibrary::pretty_print( std::string atom1, std::string atom2 ) const {
	mm_bondlength_library_citer_pair temppair = this->lookup(atom1, atom2);
	for ( mm_bondlength_library_citer i = temppair.first, e = temppair.second; i != e; ++i ) {
		TR << (i->first).key1() << "\t"
			<< (i->first).key2() << "\t"
			<< (i->second).key1() << "\t"
			<< (i->second).key2() << "\t"
			<< std::endl;
	}
}


} // namespace mm
} // namespace scoring
} // namespace core
