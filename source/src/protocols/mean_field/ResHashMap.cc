// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    ResHashMap.cc

/// @brief   Method definitions for ResHashMap.
/// @author  arubenstein

// Unit headers
#include <protocols/mean_field/ResHashMap.hh>

// Package headers


// Project headers

// Numeric headers


// C++ headers
#include <iostream>


namespace protocols {
namespace mean_field {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////
// Default constructor
/// @details
ResHashMap::ResHashMap() : utility::pointer::ReferenceCount(),
	hash_( Size( 1000 ) ),
	last_ind_assigned_( Size ( 0 ) )
{}

// Destructor
ResHashMap::~ResHashMap() {}


// Standard Rosetta methods ////////////////////////////////////////////////////
// General methods
void
ResHashMap::show( std::ostream & output ) const
{
	output << "ResHashMap of size " << hash_.size();
}

/// @details Checks if residue res is in hash map, if it is returns its RotamerIndex
/// @details If it is not in the hash map, inserts the residue res with a RotamerIndex of ind and returns ind
protocols::mean_field::ResHashMap::RotamerIndex
ResHashMap::attempt_insert ( conformation::ResidueCOP res )
{
	RotamerIndex rot_ind = get_rot_ind( res );
	if ( rot_ind == -1 ) {
		rot_ind = ++last_ind_assigned_;
		hash_.insert( std::make_pair( res , rot_ind ) );
	}
	return rot_ind;
}

/// @details if res is already in the hash map, returns its rotamer index
/// @details otherwise, returns -1
protocols::mean_field::ResHashMap::RotamerIndex
ResHashMap::get_rot_ind ( conformation::ResidueCOP res ) const
{
	RotamerIndex rot_ind = Size ( -1 );

	ResHash::const_iterator result = hash_.find( res );

	if ( result != hash_.end() ) {
		rot_ind = result->second;
	}

	return rot_ind;

}

// Friend methods //////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that ResHashMap can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, ResHashMap const & object_to_output )
{
	object_to_output.show(output);
	return output;
}

}  // namespace mean_field
}  // namespace protocols
