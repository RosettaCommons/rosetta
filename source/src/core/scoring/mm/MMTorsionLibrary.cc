// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMTorsionLibary.cc
/// @brief  Molecular mechanics torsion library
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers
#include <core/scoring/mm/MMTorsionLibrary.hh>

// Project headers
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/keys/Key4Tuple.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers
#include <numeric/conversions.hh>

// External headers
#include <boost/lexical_cast.hpp>

// C++ headers
#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include <utility/assert.hh>


namespace core {
namespace scoring {
namespace mm {

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "core.mm.MMTorsionLibrary" );

MMTorsionLibrary::~MMTorsionLibrary() {}

/// @details Constructs a MMTorsionLibrary instance from a filename string and constant access pointer to an MMAtomTypeSet
MMTorsionLibrary::MMTorsionLibrary(
	std::string filename,
	core::chemical::MMAtomTypeSetCOP mm_atom_set
)
{
	mm_atom_set_ = core::chemical::MMAtomTypeSetCAP( mm_atom_set );

	// read the file
	std::string line;
	utility::vector1< std::string > lines;
	std::ifstream data( filename.c_str() );

	while ( getline( data, line ) ) {
		if ( line.size() < 1 || line[0] == '#' ) continue; // comment or blank lines
		lines.push_back( line );
	}

	// add the torsion params
	for ( Size i = 1; i <= lines.size(); ++i ) {

		std::istringstream l( lines[i] );

		// get four atom type strings
		std::string atom_type_string_1, atom_type_string_2,
			atom_type_string_3, atom_type_string_4;
		l >> atom_type_string_1 >> atom_type_string_2
			>> atom_type_string_3 >> atom_type_string_4;

		// get atom-type_index from atom set
		int atom_type_int1 = mm_atom_set->atom_type_index( atom_type_string_1 );
		int atom_type_int2 = mm_atom_set->atom_type_index( atom_type_string_2 );
		int atom_type_int3 = mm_atom_set->atom_type_index( atom_type_string_3 );
		int atom_type_int4 = mm_atom_set->atom_type_index( atom_type_string_4 );

		// get k_theta, multiplicty, and minimum
		Real k_theta, minimum; int multiplicity;
		l >> k_theta >>  multiplicity >> minimum;
		minimum = numeric::conversions::radians( minimum );

		// add to correct library
		if ( atom_type_string_1 == "X" || atom_type_string_2 == "X" || atom_type_string_3 == "X" || atom_type_string_4 == "X" ) {
			wildcard_mm_torsion_library_.insert( std::make_pair(
				mm_torsion_atom_quad( atom_type_int1, atom_type_int2, atom_type_int3, atom_type_int4 ),
				mm_torsion_param_set( k_theta, multiplicity, minimum ) ) );
		} else {
			fully_assigned_mm_torsion_library_.insert( std::make_pair(
				mm_torsion_atom_quad( atom_type_int1, atom_type_int2, atom_type_int3, atom_type_int4 ),
				mm_torsion_param_set( k_theta, multiplicity, minimum ) ) );
		}
	}

	/// apl -- add "no-op" pair for virtual atoms
	int const virt_type = mm_atom_set->atom_type_index("VIRT");
	Real const no_op_k_theta( 0.0 );
	Real const no_op_minimum( 0.0 );
	int const no_op_multiplicity( 0 );
	fully_assigned_mm_torsion_library_.insert( std::make_pair(
		mm_torsion_atom_quad( virt_type, virt_type, virt_type, virt_type ),
		mm_torsion_param_set( no_op_k_theta, no_op_multiplicity, no_op_minimum ) ));

	// print number torsion params added
	TR << "MM torsion sets added fully assigned: " << fully_assigned_mm_torsion_library_.size()
		<< "; wildcard: " << wildcard_mm_torsion_library_.size() << " and 1 virtual parameter."
		<< std::endl;
}

/// @details The lookup function returns a pair of iterators to the first and last element in the multimap library that
/// corespond to the set(s) of mm params that corespond to the 4 mm atom type indices. If non are found it exits.
mm_torsion_library_citer_pair
MMTorsionLibrary::lookup(
	int atom1,
	int atom2,
	int atom3,
	int atom4) const
{
	static std::string const x_string = "X";
	static std::string const virt_string = "VIRT";

	// forward
	if ( fully_assigned_mm_torsion_library_.count(
			mm_torsion_atom_quad( atom1, atom2, atom3, atom4 ) ) ) {
		return fully_assigned_mm_torsion_library_.equal_range(
			mm_torsion_atom_quad( atom1, atom2, atom3, atom4 ) );
	} else if ( fully_assigned_mm_torsion_library_.count(
			mm_torsion_atom_quad( atom4, atom3, atom2, atom1 ) ) ) {
		// backward
		return fully_assigned_mm_torsion_library_.equal_range(
			mm_torsion_atom_quad( atom4, atom3, atom2, atom1 ) );
	}

	core::chemical::MMAtomTypeSetCOP mm_atom_set( mm_atom_set_ );
	int const virt_atom_type = mm_atom_set->atom_type_index( virt_string );

	// Virtual atoms get no mm-parameters.  Return the no-op torsion object
	if ( atom1 == virt_atom_type ||
			atom2 == virt_atom_type ||
			atom3 == virt_atom_type ||
			atom4 == virt_atom_type ) {
		return fully_assigned_mm_torsion_library_.equal_range(
			mm_torsion_atom_quad( virt_atom_type, virt_atom_type, virt_atom_type, virt_atom_type ));
	}


	int const wild_atom_type = mm_atom_set->atom_type_index( x_string );

	// check to see if there are any wildcard params want to use the wild card with the most atoms defined
	// this is totally ugly and needs to be changed, maybe some fast hash table based lookup would be faster
	// all of this looking through maps

	// 1 wild card forward (XAAA, AXAA, AAXA, AAAX)
	if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, atom2, atom3, atom4 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, atom2, atom3, atom4 ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom1, wild_atom_type, atom3, atom4 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom1, wild_atom_type, atom3, atom4 ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom1, atom2, wild_atom_type, atom4 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom1, atom2, wild_atom_type, atom4 ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom1, atom2, atom3, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom1, atom2, atom3, wild_atom_type ) ); }
	// 1 wild card backward  (XAAA, AXAA, AAXA, AAAX)
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, atom3, atom2, atom1 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, atom3, atom2, atom1 ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom4, wild_atom_type, atom2, atom1 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom4, wild_atom_type, atom2, atom1 ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom4, atom3, wild_atom_type, atom1 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom4, atom3, wild_atom_type, atom1 ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom4, atom3, atom2, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom4, atom3, atom2, wild_atom_type ) ); }

	// 2 wild card forward (AAXX, AXXA, XXAA, AXAX, XAXA, XAAX)
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom1, atom2, wild_atom_type, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom1, atom2, wild_atom_type, wild_atom_type ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom1, wild_atom_type, wild_atom_type, atom4 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom1, wild_atom_type, wild_atom_type, atom4 ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, wild_atom_type, atom3, atom4 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, wild_atom_type, atom3, atom4 ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom1, wild_atom_type, atom3, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom1, wild_atom_type, atom3, wild_atom_type ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, atom2, wild_atom_type, atom4 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, atom2, wild_atom_type, atom4 ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, atom2, atom3, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, atom2, atom3, wild_atom_type ) ); }
	// 2 wild card backward (AAXX, AXXA, XXAA, AXAX, XAXA, XAAX)
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom4, atom3, wild_atom_type, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom4, atom3, wild_atom_type, wild_atom_type ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom4, wild_atom_type, wild_atom_type, atom1 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom4, wild_atom_type, wild_atom_type, atom1 ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, wild_atom_type, atom2, atom1 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, wild_atom_type, atom2, atom1 ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom4, wild_atom_type, atom3, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom4, wild_atom_type, atom3, wild_atom_type ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, atom3, wild_atom_type, atom1 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, atom3, wild_atom_type, atom1 ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, atom3, atom2, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, atom3, atom2, wild_atom_type ) ); }

	// 3 wildcard forward (AXXX, XAXX, XXAX, XXXA)
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom1, wild_atom_type, wild_atom_type, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom1, wild_atom_type, wild_atom_type, wild_atom_type ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, atom2, wild_atom_type, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, atom2, wild_atom_type, wild_atom_type ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, wild_atom_type, atom3, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, wild_atom_type, atom3, wild_atom_type ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, wild_atom_type, wild_atom_type, atom4 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, wild_atom_type, wild_atom_type, atom4 ) ); }
	// 3 wildcard backward (AXXX, XAXX, XXAX, XXXA)
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( atom4, wild_atom_type, wild_atom_type, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( atom4, wild_atom_type, wild_atom_type, wild_atom_type ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, atom3, wild_atom_type, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, atom3, wild_atom_type, wild_atom_type ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, wild_atom_type, atom2, wild_atom_type ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, wild_atom_type, atom2, wild_atom_type ) ); }
	else if ( wildcard_mm_torsion_library_.count( mm_torsion_atom_quad( wild_atom_type, wild_atom_type, wild_atom_type, atom1 ) ) ) { return wildcard_mm_torsion_library_.equal_range( mm_torsion_atom_quad( wild_atom_type, wild_atom_type, wild_atom_type, atom1 ) ); }

	TR << "No parameters for "
		<< (*mm_atom_set)[atom1].name() << "-"
		<< (*mm_atom_set)[atom2].name() << "-"
		<< (*mm_atom_set)[atom3].name() << "-"
		<< (*mm_atom_set)[atom4].name() << std::endl;

	/// Arriving here, we've failed to find essential parameters
	utility_exit_with_message("COULD NOT FIND TORSION PARAMS FOR " +
		boost::lexical_cast<std::string>(atom1) + " " +
		boost::lexical_cast<std::string>(atom2) + " " +
		boost::lexical_cast<std::string>(atom3) + " " +
		boost::lexical_cast<std::string>(atom4) );
	return mm_torsion_library_citer_pair();  //< meaningless, just for removing gcc warning.
}

/// @details
mm_torsion_library_citer_pair
MMTorsionLibrary::lookup
(
	std::string atom1,
	std::string atom2,
	std::string atom3,
	std::string atom4
) const
{
	core::chemical::MMAtomTypeSetCOP mm_atom_set( mm_atom_set_ );
	return (*this).lookup( mm_atom_set->atom_type_index( atom1 ),
		mm_atom_set->atom_type_index( atom2 ),
		mm_atom_set->atom_type_index( atom3 ),
		mm_atom_set->atom_type_index( atom4 ) );
}

} // namespace mm
} // namespace scoring
} // namespace core
