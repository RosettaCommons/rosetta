// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/CrystInfo.cc
/// @brief  Auto-generated serialization template functions
/// @author Andrew Leaver-Fay (aleavefay@gmail.com)

// Unit headers
#include <core/io/CrystInfo.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>

/// @brief Automatically generated serialization method
template< class Archive >
void
core::io::CrystInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( A_ ) ); // Real
	arc( CEREAL_NVP( B_ ) ); // Real
	arc( CEREAL_NVP( C_ ) ); // Real
	arc( CEREAL_NVP( alpha_ ) ); // Real
	arc( CEREAL_NVP( beta_ ) ); // Real
	arc( CEREAL_NVP( gamma_ ) ); // Real
	arc( CEREAL_NVP( spacegroup_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::io::CrystInfo::load( Archive & arc ) {
	arc( A_ ); // Real
	arc( B_ ); // Real
	arc( C_ ); // Real
	arc( alpha_ ); // Real
	arc( beta_ ); // Real
	arc( gamma_ ); // Real
	arc( spacegroup_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::io::CrystInfo );
#endif // SERIALIZATION

