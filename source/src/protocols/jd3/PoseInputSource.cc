// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/PoseInputSource.cc
/// @brief  Definition of the %PoseInputSource class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

//unit headers
#include <protocols/jd3/PoseInputSource.hh>

//project headers
#include <core/pose/Pose.fwd.hh>
#include <utility>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {

PoseInputSource::PoseInputSource() :
	origin_( "unknown" )
{}

PoseInputSource::PoseInputSource( std::string  origin ) :
	origin_(std::move( origin ))
{}
PoseInputSource::~PoseInputSource() = default;

bool PoseInputSource::operator == ( PoseInputSource const & rhs ) const
{
	return
		origin_ == rhs.origin_ &&
		input_tag_ == rhs.input_tag_ &&
		string_string_map_ == rhs.string_string_map_;
}

bool PoseInputSource::operator != ( PoseInputSource const & rhs ) const
{
	return ! ( *this == rhs );
}

bool PoseInputSource::operator < ( PoseInputSource const & rhs ) const
{
	if ( origin_ <  rhs.origin_ ) return true;
	if ( origin_ != rhs.origin_ ) return false;
	if ( input_tag_ <  rhs.input_tag_ ) return true;
	if ( input_tag_ != rhs.input_tag_ ) return false;
	if ( string_string_map_ < rhs.string_string_map_ ) return true;
	return false;
}

std::string const & PoseInputSource::input_tag() const { return input_tag_; }
PoseInputSource::StringStringMap const &
PoseInputSource::string_string_map() const { return string_string_map_; }
std::string const & PoseInputSource::origin() const { return origin_; }

void PoseInputSource::input_tag( std::string const & setting ) { input_tag_ = setting; }
void PoseInputSource::store_string_pair( std::string const & key, std::string const & value ) { string_string_map_[ key ] = value; }
void PoseInputSource::origin( std::string const & setting ) { origin_ = setting; }

} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::PoseInputSource::save( Archive & arc ) const {
	arc( CEREAL_NVP( origin_ ) ); // std::string
	arc( CEREAL_NVP( input_tag_ ) ); // std::string
	arc( CEREAL_NVP( string_string_map_ ) ); // StringStringMap
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::PoseInputSource::load( Archive & arc ) {
	arc( origin_ ); // std::string
	arc( input_tag_ ); // std::string
	arc( string_string_map_ ); // StringStringMap
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::PoseInputSource );
CEREAL_REGISTER_TYPE( protocols::jd3::PoseInputSource )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_PoseInputSource )
#endif // SERIALIZATION
