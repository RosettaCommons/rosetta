// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.

/// @file   utility/options/Option.cc
/// @brief  Auto-generated serialization template functions
/// @author Andrew Leaver-Fay (aleavefay@gmail.com)

// Unit headers
#include <utility/options/Option.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>

/// @brief Automatically generated serialization method
template< class Archive >
void
utility::options::Option::save( Archive & arc ) const {
	arc( CEREAL_NVP( is_group_ ) ); // _Bool
#ifdef MULTI_THREADED
	bool been_accessed( been_accessed_.load() );
	arc( CEREAL_NVP( been_accessed ) ); // _Bool
#else
	arc( CEREAL_NVP( been_accessed_ ) ); // _Bool
#endif
	arc( CEREAL_NVP( restricted_access_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
utility::options::Option::load( Archive & arc ) {
	arc( is_group_ ); // _Bool
	bool been_accessed;
	arc( been_accessed ); // _Bool
#ifdef MULTI_THREADED
	been_accessed_.store(been_accessed);
#else
	been_accessed_ = been_accessed;
#endif
	arc( restricted_access_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( utility::options::Option );
CEREAL_REGISTER_TYPE( utility::options::Option )

CEREAL_REGISTER_DYNAMIC_INIT( utility_options_Option )

#endif // SERIALIZATION
