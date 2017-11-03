// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/DOF_ID.cc
/// @brief  Kinematics DOF identifier class
/// @author Phil Bradley


// Unit headers
#include <core/id/DOF_ID.hh>

// C++ headers

#include <iostream>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

#endif // SERIALIZATION


namespace core {
namespace id {


/// @brief stream << DOF_ID
std::ostream &
operator <<( std::ostream & os, DOF_ID const & a )
{
	os << " atom_id= " << a.atom_id() << " type= " << a.type() << ' ';
	return os;
}


} // namespace id
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::id::DOF_ID::save( Archive & arc ) const {
	arc( CEREAL_NVP( atom_id_ ) ); // class core::id::AtomID
	arc( CEREAL_NVP( type_ ) ); // enum core::id::DOF_Type
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::id::DOF_ID::load( Archive & arc ) {
	arc( atom_id_ ); // class core::id::AtomID
	arc( type_ ); // enum core::id::DOF_Type
}

SAVE_AND_LOAD_SERIALIZABLE( core::id::DOF_ID );
#endif // SERIALIZATION

