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
/// @author Modified by Vikram K. Mulligan (vmulligan@flatironinstitute.org) for thread-safety.

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

#ifdef MULTI_THREADED
	bool is_group( is_group_ );
	bool been_accessed( been_accessed_ );
	bool restricted_access( restricted_access_ );
	arc( CEREAL_NVP( is_group ) ); // _Bool
	arc( CEREAL_NVP( been_accessed ) ); // _Bool
	arc( CEREAL_NVP( restricted_access ) ); // _Bool
#else
	arc( CEREAL_NVP( is_group_ ) ); // _Bool
	arc( CEREAL_NVP( been_accessed_ ) ); // _Bool
	arc( CEREAL_NVP( restricted_access_ ) ); // _Bool
#endif
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
utility::options::Option::load( Archive & arc ) {

#ifdef MULTI_THREADED
	bool is_group, been_accessed, restricted_access;
	arc( is_group ); // _Bool
	arc( been_accessed ); // _Bool
	arc( restricted_access ); // _Bool
	is_group_ = is_group;
	been_accessed_ = been_accessed;
	restricted_access_ = restricted_access;
#else
	arc( is_group_ ); // _Bool
	arc( been_accessed_);
	arc( restricted_access_ ); // _Bool
#endif
}

SAVE_AND_LOAD_SERIALIZABLE( utility::options::Option );
CEREAL_REGISTER_TYPE( utility::options::Option )

CEREAL_REGISTER_DYNAMIC_INIT( utility_options_Option )

#endif // SERIALIZATION
