// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.

/// @file   core/chemical/ResidueConnection.cc
/// @brief  Auto-generated serialization template functions
/// @author Andrew Leaver-Fay (aleavefay@gmail.com)

// Unit headers
#include <core/chemical/ResidueConnection.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#endif // SERIALIZATION

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ResidueConnection::save( Archive & arc ) const {
	arc( CEREAL_NVP( atomno_ ) ); // int
	arc( CEREAL_NVP( icoor_ ) ); // class core::chemical::AtomICoor
	arc( CEREAL_NVP( index_ ) ); // int
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ResidueConnection::load( Archive & arc ) {
	arc( atomno_ ); // int
	arc( icoor_ ); // class core::chemical::AtomICoor
	arc( index_ ); // int
}
SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ResidueConnection );
#endif // SERIALIZATION
