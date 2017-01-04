// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/deallocation/InputPoseDeallocationMessage.cc
/// @brief  The implementation for class protocols::jd3::deallocation::InputPoseDeallocationMessage's methods
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/deallocation/InputPoseDeallocationMessage.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace deallocation {

InputPoseDeallocationMessage::InputPoseDeallocationMessage() :
	DeallocationMessage( input_pose_deallocation_msg ),
	pose_id_( 0 )
{}

InputPoseDeallocationMessage::InputPoseDeallocationMessage( core::Size pose_id ) :
	DeallocationMessage( input_pose_deallocation_msg ),
	pose_id_( pose_id )
{}

InputPoseDeallocationMessage::~InputPoseDeallocationMessage() {}

core::Size InputPoseDeallocationMessage::pose_id() const { return pose_id_; }

void InputPoseDeallocationMessage::pose_id( core::Size setting )
{
	pose_id_ = setting;
}

} // namespace deallocation
} // namespace jd3
} // namespace protocols


#ifdef    SERIALIZATION

template< class Archive >
void
protocols::jd3::deallocation::InputPoseDeallocationMessage::save( Archive & arc ) const {
	arc( cereal::base_class< DeallocationMessage >( this ) );
	arc( CEREAL_NVP( pose_id_ ) );
}

template< class Archive >
void
protocols::jd3::deallocation::InputPoseDeallocationMessage::load( Archive & arc ) {
	arc( cereal::base_class< DeallocationMessage >( this ) );
	arc( pose_id_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::deallocation::InputPoseDeallocationMessage );
CEREAL_REGISTER_TYPE( protocols::jd3::deallocation::InputPoseDeallocationMessage )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_deallocation_InputPoseDeallocationMessage )
#endif // SERIALIZATION
