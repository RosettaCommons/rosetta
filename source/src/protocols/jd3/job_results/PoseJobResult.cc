// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/job_results/PoseJobResult.cc
/// @brief  The class method definitions for PoseJobResult
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/job_results/PoseJobResult.hh>

// Project headers
#include <core/pose/Pose.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace job_results {

PoseJobResult::PoseJobResult() = default;

PoseJobResult::PoseJobResult( core::pose::PoseOP setting ) :
	pose_( std::move( setting ) )
{}

PoseJobResult::~PoseJobResult() = default;

JobStatus
PoseJobResult::status() const{
	return jd3_job_status_success;
}

core::pose::PoseOP
PoseJobResult::pose() {
	return pose_;
}

core::pose::PoseCOP
PoseJobResult::pose() const {
	return pose_;
}

void
PoseJobResult::pose( core::pose::PoseOP setting ) {
	pose_ = setting;
}


} // namespace job_results
} // namespace jd3
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::job_results::PoseJobResult::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::JobResult >( this ) );
	arc( CEREAL_NVP( pose_ ) ); // core::pose::PoseOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::job_results::PoseJobResult::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::JobResult >( this ) );
	arc( pose_ ); // core::pose::PoseOP
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::job_results::PoseJobResult );
CEREAL_REGISTER_TYPE( protocols::jd3::job_results::PoseJobResult )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_job_results_PoseJobResult )
#endif // SERIALIZATION
