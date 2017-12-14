// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_outputters/PoseOutputSpecification.cc
/// @brief  Method definitions of the %PoseOutputSpecification class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/pose_outputters/PoseOutputSpecification.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace pose_outputters {

PoseOutputSpecification::PoseOutputSpecification() = default;

PoseOutputSpecification::PoseOutputSpecification(
	JobResultID const & result_id,
	JobOutputIndex const & output_index
) :
	OutputSpecification( result_id, output_index )
{}

PoseOutputSpecification::~PoseOutputSpecification() = default;

std::string const & PoseOutputSpecification::outputter_type() const
{
	return outputter_type_;
}

void PoseOutputSpecification::outputter_type( std::string const & setting )
{
	outputter_type_ = setting;
}

} // namespace pose_outputters
} // namespace jd3
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::pose_outputters::PoseOutputSpecification::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::output::OutputSpecification >( this ) );
	arc( CEREAL_NVP( outputter_type_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::pose_outputters::PoseOutputSpecification::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::output::OutputSpecification >( this ) );
	arc( outputter_type_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::pose_outputters::PoseOutputSpecification );
CEREAL_REGISTER_TYPE( protocols::jd3::pose_outputters::PoseOutputSpecification )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_pose_outputters_PoseOutputSpecification )
#endif // SERIALIZATION
