// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/pdb//Remarks.cc
/// @brief  Auto-generated serialization template functions
/// @author Andrew Leaver-Fay (aleavefay@gmail.com)

// Unit headers
#include <core/io/Remarks.hh>

#ifdef    SERIALIZATION

// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>


/// @brief Automatically generated serialization method
template< class Archive >
void
core::io::RemarkInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( num ) ); // int
	arc( CEREAL_NVP( value ) ); // String
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::io::RemarkInfo::load( Archive & arc ) {
	arc( num ); // int
	arc( value ); // String
}

SAVE_AND_LOAD_SERIALIZABLE( core::io::RemarkInfo );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::io::Remarks::save( Archive & arc ) const {
	arc( cereal::base_class< std::vector< core::io::RemarkInfo > >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::io::Remarks::load( Archive & arc ) {
	arc( cereal::base_class< std::vector< core::io::RemarkInfo > >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::io::Remarks );
CEREAL_REGISTER_TYPE( core::io::Remarks )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_Remarks )
#endif // SERIALIZATION

