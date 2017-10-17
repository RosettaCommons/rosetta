// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/JobOutputIndex.cc
/// @brief  Auto-generated serialization template functions
/// @author Andrew Leaver-Fay (aleavefay@gmail.com)

#include <protocols/jd3/JobOutputIndex.hh>


// Unit headers
#include <protocols/jd3/JobOutputIndex.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif

// Empty body of this .cc except for serialization headers

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::JobOutputIndex::save( Archive & arc ) const {
	arc( CEREAL_NVP( primary_output_index ) ); // core::Size
	arc( CEREAL_NVP( n_primary_outputs ) ); // core::Size
	arc( CEREAL_NVP( secondary_output_index ) ); // core::Size
	arc( CEREAL_NVP( n_secondary_outputs ) ); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::JobOutputIndex::load( Archive & arc ) {
	arc( primary_output_index ); // core::Size
	arc( n_primary_outputs ); // core::Size
	arc( secondary_output_index ); // core::Size
	arc( n_secondary_outputs ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::JobOutputIndex );
CEREAL_REGISTER_TYPE( protocols::jd3::JobOutputIndex )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_JobOutputIndex )
#endif // SERIALIZATION
