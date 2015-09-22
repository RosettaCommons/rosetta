// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMBondAngleLibary.cc
/// @brief  Molecular mechanics bond angle library
/// @author Colin A. Smith (colin.smith@ucsf.edu)

// Unit headers
#include <core/scoring/mm/MMBondAngleLibrary.hh>

// Project headers
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key2Tuple.hh>
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
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility/assert.hh>


namespace core {
namespace scoring {
namespace mm {

/// @details Auto-generated virtual destructor
MMBondAngleLibrary::~MMBondAngleLibrary() {}

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "core.mm.MMBondAngleLibrary" );

/// @details Construct a MMBondAngleLibrary instant from a filename string and constant access pointer to an MMAtomTypeSet
MMBondAngleLibrary::MMBondAngleLibrary(
	std::string filename,
	core::chemical::MMAtomTypeSetCAP mm_atom_set_ap
)
{
	mm_atom_set_ = mm_atom_set_ap;

	core::chemical::MMAtomTypeSetCOP mm_atom_set( mm_atom_set_ );

	// read the file
	std::string line;
	utility::vector1< std::string > lines;
	std::ifstream data( filename.c_str() );

	bool in_bonds_section = false;
	while ( getline( data, line ) ) {
		if ( line.size() < 1 || line[0] == '!' || line[0] == ' ' ) continue; // comment or blank lines
		if ( line == "DIHEDRALS" ) in_bonds_section = false;
		if ( in_bonds_section ) lines.push_back( line );
		if ( line == "ANGLES" ) in_bonds_section = true;
	}

	// add the torsion params
	for ( Size i = 1; i <= lines.size(); ++i ) {

		std::istringstream l( lines[i] );

		// get four atom type strings
		std::string atom_type_string_1, atom_type_string_2,
			atom_type_string_3;
		l >> atom_type_string_1 >> atom_type_string_2
			>> atom_type_string_3;

		// skip the parameters if any of the atom types don't exist
		if ( ! mm_atom_set->contains_atom_type( atom_type_string_1 ) ) continue;
		if ( ! mm_atom_set->contains_atom_type( atom_type_string_2 ) ) continue;
		if ( ! mm_atom_set->contains_atom_type( atom_type_string_3 ) ) continue;

		// get atom-type_index from atom set
		int atom_type_int1 = mm_atom_set->atom_type_index( atom_type_string_1 );
		int atom_type_int2 = mm_atom_set->atom_type_index( atom_type_string_2 );
		int atom_type_int3 = mm_atom_set->atom_type_index( atom_type_string_3 );

		// get k_theta and minimum
		Real k_theta, minimum;
		l >> k_theta >> minimum;
		minimum = numeric::conversions::radians( minimum );

		//TR << atom_type_string_1 << "\t" << atom_type_string_2 << "\t" << atom_type_string_2 << "\t" << k_theta << "\t"
		//   << minimum << std::endl;

		// add to correct library
		if ( atom_type_string_1 == "X" && atom_type_string_3 == "X" ) {
			wildcard_mm_bondangle_library_.insert( std::make_pair(
				mm_bondangle_atom_tri( atom_type_int1, atom_type_int2, atom_type_int3 ),
				mm_bondangle_param_set( k_theta, minimum ) ) );
		} else {
			fully_assigned_mm_bondangle_library_.insert( std::make_pair(
				mm_bondangle_atom_tri( atom_type_int1, atom_type_int2, atom_type_int3 ),
				mm_bondangle_param_set( k_theta, minimum ) ) );
		}
	}

	/// apl -- add "no-op" pair for virtual atoms
	int const virt_type = mm_atom_set->atom_type_index("VIRT");
	Real const no_op_k_theta( 0.0 );
	Real const no_op_minimum( 0.0 );
	fully_assigned_mm_bondangle_library_.insert( std::make_pair(
		mm_bondangle_atom_tri( virt_type, virt_type, virt_type ),
		mm_bondangle_param_set( no_op_k_theta, no_op_minimum ) ));


	// print number torsion params added
	TR << "MM bond angle sets added fully assigned: " << fully_assigned_mm_bondangle_library_.size()
		<< "; wildcard: " << wildcard_mm_bondangle_library_.size() << " and 1 virtual parameter."
		<< std::endl;
}

mm_bondangle_library_citer_pair
MMBondAngleLibrary::lookup (
	int atom1,
	int atom2,
	int atom3) const
{
	static std::string const x_string = "X";
	static std::string const virt_string = "VIRT";

	if ( fully_assigned_mm_bondangle_library_.count(
			mm_bondangle_atom_tri( atom1, atom2, atom3 ) ) ) {
		// forward
		return fully_assigned_mm_bondangle_library_.equal_range(
			mm_bondangle_atom_tri( atom1, atom2, atom3 ) );
	} else if ( fully_assigned_mm_bondangle_library_.count(
			mm_bondangle_atom_tri( atom3, atom2, atom1 ) ) ) {
		// backward
		return fully_assigned_mm_bondangle_library_.equal_range(
			mm_bondangle_atom_tri( atom3, atom2, atom1 ) );
	}

	core::chemical::MMAtomTypeSetCOP mm_atom_set( mm_atom_set_ );
	int const virt_atom_type = mm_atom_set->atom_type_index( virt_string );

	// Virtual atoms get no mm-parameters.  Return the no-op torsion object
	if ( atom1 == virt_atom_type ||
			atom2 == virt_atom_type ||
			atom3 == virt_atom_type ) {
		return fully_assigned_mm_bondangle_library_.equal_range(
			mm_bondangle_atom_tri( virt_atom_type, virt_atom_type, virt_atom_type ));
	}


	int const wild_atom_type = mm_atom_set->atom_type_index( x_string );

	if ( wildcard_mm_bondangle_library_.count(
			mm_bondangle_atom_tri( wild_atom_type, atom2, wild_atom_type ) ) ) {
		// wildcard 1 & 3
		return wildcard_mm_bondangle_library_.equal_range(
			mm_bondangle_atom_tri( wild_atom_type, atom2, wild_atom_type ) );
	}

	TR << "No parameters for " << (*mm_atom_set)[atom1].name() << "-" << (*mm_atom_set)[atom2].name() << "-"
		<< (*mm_atom_set)[atom3].name() << std::endl;
	utility_exit_with_message("COULD NOT FIND BOND ANGLE PARAMS FOR " +
		boost::lexical_cast<std::string>(atom1) + " " +
		boost::lexical_cast<std::string>(atom2) + " " +
		boost::lexical_cast<std::string>(atom3) );

	return mm_bondangle_library_citer_pair();  //< meaningless, just for removing gcc warning.
}

mm_bondangle_library_citer_pair
MMBondAngleLibrary::lookup
(
	std::string atom1,
	std::string atom2,
	std::string atom3
) const
{
	core::chemical::MMAtomTypeSetCOP mm_atom_set( mm_atom_set_ );
	return (*this).lookup( mm_atom_set->atom_type_index( atom1 ),
		mm_atom_set->atom_type_index( atom2 ),
		mm_atom_set->atom_type_index( atom3 ) );
}

void
MMBondAngleLibrary::pretty_print() const
{
	// for each key print out its value
	for ( mm_bondangle_library_citer i = fully_assigned_mm_bondangle_library_.begin(),
			e = fully_assigned_mm_bondangle_library_.end(); i != e; ++i ) {
		TR << (i->first).key1() << "\t"
			<< (i->first).key2() << "\t"
			<< (i->first).key3() << "\t"
			<< (i->second).key1() << "\t"
			<< (i->second).key2() << "\t"
			<< std::endl;
	}

	for ( mm_bondangle_library_citer i = wildcard_mm_bondangle_library_.begin(),
			e = wildcard_mm_bondangle_library_.end(); i != e; ++i ) {
		TR << (i->first).key1() << "\t"
			<< (i->first).key2() << "\t"
			<< (i->first).key3() << "\t"
			<< (i->second).key1() << "\t"
			<< (i->second).key2() << "\t"
			<< std::endl;
	}
}

void
MMBondAngleLibrary::pretty_print(  int atom1, int atom2, int atom3 ) const
{
	mm_bondangle_library_citer_pair temppair = (*this).lookup(atom1, atom2, atom3);
	for ( mm_bondangle_library_citer i = temppair.first, e = temppair.second; i != e; ++i ) {
		TR << (i->first).key1() << "\t"
			<< (i->first).key2() << "\t"
			<< (i->first).key3() << "\t"
			<< (i->second).key1() << "\t"
			<< (i->second).key2() << "\t"
			<< std::endl;
	}
}

void
MMBondAngleLibrary::pretty_print( std::string atom1, std::string atom2, std::string atom3 ) const
{
	mm_bondangle_library_citer_pair temppair = (*this).lookup(atom1, atom2, atom3);
	for ( mm_bondangle_library_citer i = temppair.first, e = temppair.second; i != e; ++i ) {
		TR << (i->first).key1() << "\t"
			<< (i->first).key2() << "\t"
			<< (i->first).key3() << "\t"
			<< (i->second).key1() << "\t"
			<< (i->second).key2() << "\t"
			<< std::endl;
	}
}

} // namespace mm
} // namespace scoring
} // namespace core
