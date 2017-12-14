// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_outputters/PDBPoseOutputSpecification.cc
/// @brief  Method definitions of the %PDBPoseOutputSpecification class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/pose_outputters/PDBPoseOutputSpecification.hh>

// Package headers
#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>

// Project headers
#include <core/io/StructFileRepOptions.hh>

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

PDBPoseOutputSpecification::PDBPoseOutputSpecification()
{
	outputter_type( PDBPoseOutputter::keyname() );
}

PDBPoseOutputSpecification::PDBPoseOutputSpecification(
	JobResultID const & result_id,
	JobOutputIndex const & output_index
) :
	PoseOutputSpecification( result_id, output_index )
{
	outputter_type( PDBPoseOutputter::keyname() );
}

PDBPoseOutputSpecification::~PDBPoseOutputSpecification() = default;


core::io::StructFileRepOptionsCOP
PDBPoseOutputSpecification::sfr_opts() const
{
	debug_assert( sfr_opts_ );
	return sfr_opts_;
}

void
PDBPoseOutputSpecification::sfr_opts( core::io::StructFileRepOptions const & setting )
{
	sfr_opts_.reset( new core::io::StructFileRepOptions( setting ) );
}

std::string
PDBPoseOutputSpecification::out_fname() const
{
	return out_fname_;
}

void
PDBPoseOutputSpecification::out_fname( std::string const & setting )
{
	out_fname_ = setting;
}



} // namespace pose_outputters
} // namespace jd3
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::pose_outputters::PDBPoseOutputSpecification::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::pose_outputters::PoseOutputSpecification >( this ) );
	arc( CEREAL_NVP( sfr_opts_ ) ); // core::io::StructFileRepOptionsOP
	arc( CEREAL_NVP( out_fname_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::pose_outputters::PDBPoseOutputSpecification::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::pose_outputters::PoseOutputSpecification >( this ) );
	arc( sfr_opts_ ); // core::io::StructFileRepOptionsOP
	arc( out_fname_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::pose_outputters::PDBPoseOutputSpecification );
CEREAL_REGISTER_TYPE( protocols::jd3::pose_outputters::PDBPoseOutputSpecification )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_pose_outputters_PDBPoseOutputSpecification )
#endif // SERIALIZATION
