// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/deallocation/ResourceDeallocationMessage.cc
/// @brief  The implementation for class protocols::jd3::deallocation::ResourceDeallocationMessage's methods
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/deallocation/ResourceDeallocationMessage.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace deallocation {

ResourceDeallocationMessage::ResourceDeallocationMessage() :
	DeallocationMessage( resource_deallocation_msg )
{}

ResourceDeallocationMessage::ResourceDeallocationMessage( std::string const & resource_name ) :
	DeallocationMessage( resource_deallocation_msg ),
	resource_name_( resource_name )
{}

ResourceDeallocationMessage::~ResourceDeallocationMessage() {}

std::string const & ResourceDeallocationMessage::resource_name() const
{
	return resource_name_;
}

void ResourceDeallocationMessage::resource_name( std::string const & setting )
{
	resource_name_ = setting;
}

} // namespace deallocation
} // namespace jd3
} // namespace protocols


#ifdef    SERIALIZATION

template< class Archive >
void
protocols::jd3::deallocation::ResourceDeallocationMessage::save( Archive & arc ) const {
	arc( cereal::base_class< DeallocationMessage >( this ) );
	arc( CEREAL_NVP( resource_name_ ) );
}

template< class Archive >
void
protocols::jd3::deallocation::ResourceDeallocationMessage::load( Archive & arc ) {
	arc( cereal::base_class< DeallocationMessage >( this ) );
	arc( resource_name_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::deallocation::ResourceDeallocationMessage );
CEREAL_REGISTER_TYPE( protocols::jd3::deallocation::ResourceDeallocationMessage )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_deallocation_ResourceDeallocationMessage )
#endif // SERIALIZATION
