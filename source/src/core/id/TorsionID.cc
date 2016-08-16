// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  Kinematics Atom identifier
/// @author Phil Bradley


// Unit header
#include <core/id/TorsionID.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#endif // SERIALIZATION


namespace core {
namespace id {

/// @brief Globals
TorsionID const BOGUS_TORSION_ID( 0, BB, 0 );

} // namespace id
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::id::TorsionID::save( Archive & arc ) const {
	arc( CEREAL_NVP( rsd_ ) ); // Size
	arc( CEREAL_NVP( type_ ) ); // enum core::id::TorsionType
	arc( CEREAL_NVP( torsion_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::id::TorsionID::load( Archive & arc ) {
	arc( rsd_ ); // Size
	arc( type_ ); // enum core::id::TorsionType
	arc( torsion_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::id::TorsionID );
#endif // SERIALIZATION

