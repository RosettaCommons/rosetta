// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_outputters/SilentFilePoseOutputSpecification.cc
/// @brief  Method definitions of the %SilentFilePoseOutputSpecification class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/pose_outputters/SilentFilePoseOutputSpecification.hh>

// Package headers
#include <protocols/jd3/pose_outputters/SilentFilePoseOutputter.hh>

// Project headers
#include <core/io/silent/SilentFileOptions.hh>

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

SilentFilePoseOutputSpecification::SilentFilePoseOutputSpecification() :
	buffer_limit_( 0 )
{
	outputter_type( SilentFilePoseOutputter::keyname() );
}

SilentFilePoseOutputSpecification::SilentFilePoseOutputSpecification(
	JobResultID const & result_id,
	JobOutputIndex const & output_index
) :
	PoseOutputSpecification( result_id, output_index ),
	buffer_limit_( 0 )
{
	outputter_type( SilentFilePoseOutputter::keyname() );
}

SilentFilePoseOutputSpecification::~SilentFilePoseOutputSpecification() {}


core::io::silent::SilentFileOptions const &
SilentFilePoseOutputSpecification::sf_opts() const
{
	debug_assert( sf_opts_ );
	return *sf_opts_;
}

void
SilentFilePoseOutputSpecification::sf_opts( core::io::silent::SilentFileOptions const & setting )
{
	sf_opts_.reset( new core::io::silent::SilentFileOptions( setting ) );
}

std::string
SilentFilePoseOutputSpecification::out_fname() const
{
	return out_fname_ + suffix_from_jd_with_sep();
}

void
SilentFilePoseOutputSpecification::out_fname( std::string const & setting )
{
	out_fname_ = setting;
}

std::string const &
SilentFilePoseOutputSpecification::pose_tag() const
{
	return pose_tag_;
}

void
SilentFilePoseOutputSpecification::pose_tag( std::string const & setting )
{
	pose_tag_ = setting;
}

core::Size SilentFilePoseOutputSpecification::buffer_limit() const
{
	return buffer_limit_;
}

void SilentFilePoseOutputSpecification::buffer_limit( core::Size setting )
{
	buffer_limit_ = setting;
}


} // namespace pose_outputters
} // namespace jd3
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::pose_outputters::SilentFilePoseOutputSpecification::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::pose_outputters::PoseOutputSpecification >( this ) );
	arc( CEREAL_NVP( sf_opts_ ) ); // core::io::silent::SilentFileOptionsOP
	arc( CEREAL_NVP( out_fname_ ) ); // std::string
	arc( CEREAL_NVP( pose_tag_ ) ); // std::string
	arc( CEREAL_NVP( buffer_limit_ ) ); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::pose_outputters::SilentFilePoseOutputSpecification::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::pose_outputters::PoseOutputSpecification >( this ) );
	arc( sf_opts_ ); // core::io::silent::SilentFileOptionsOP
	arc( out_fname_ ); // std::string
	arc( pose_tag_ ); // std::string
	arc( buffer_limit_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::pose_outputters::SilentFilePoseOutputSpecification );
CEREAL_REGISTER_TYPE( protocols::jd3::pose_outputters::SilentFilePoseOutputSpecification )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_pose_outputters_SilentFilePoseOutputSpecification )
#endif // SERIALIZATION
