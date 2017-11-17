// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/integer_mapping.cc
/// @brief  A set of useful classes to map between two enumerations.  So far, only a subset mapping is implemented.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <utility/integer_mapping.hh>

// Package headers
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

// Platform headers
#include <platform/types.hh>

namespace utility {

platform::Size const subset_mapping::UNMAPPED = 0;

subset_mapping::subset_mapping() {}
subset_mapping::subset_mapping( platform::Size source_enumeration_size ) : src_2_dst_( source_enumeration_size, UNMAPPED ) {}
subset_mapping::subset_mapping( subset_mapping const & src ) : ReferenceCount(), src_2_dst_( src.src_2_dst_ ), dst_2_src_( src.dst_2_src_ ) {}
subset_mapping & subset_mapping::operator = ( subset_mapping const & rhs )
{
	if ( this != &rhs ) {
		src_2_dst_ = rhs.src_2_dst_;
		dst_2_src_ = rhs.dst_2_src_;
	}
	return *this;
}

subset_mapping::~subset_mapping() {}

void subset_mapping::set_source_size( platform::Size src_size )
{
	src_2_dst_.resize( src_size );
	std::fill( src_2_dst_.begin(), src_2_dst_.end(), UNMAPPED );
}

void subset_mapping::reserve_destination_size( platform::Size dst_size )
{
	dst_2_src_.reserve( dst_size );
}

/// @details Prevent out-of-bounds assignment and the overwriting of previously
/// mapped ids.  Throws excn::EXCN_Msg_Exceptions
void subset_mapping::set_next_correspondence( platform::Size source_id )
{
	if ( source_id > src_2_dst_.size() || source_id == 0 ) {
		throw CREATE_EXCEPTION(excn::Exception, "subset_mapping::set_next_correspondence "
			"recieved an out-of-bounds source id (" + to_string( source_id ) + ") "
			"with a source-enumeration size of " + to_string( src_2_dst_.size() ) );
	}
	if ( src_2_dst_[ source_id ] != UNMAPPED ) {
		throw CREATE_EXCEPTION(excn::Exception, "subset_mapping::set_next_correspondence "
			"recieved an already-mapped source id (" + to_string( source_id ) + ") "
			"which had been previously assigned to destination id " + to_string( src_2_dst_[ source_id ] ) );
	}
	dst_2_src_.push_back( source_id );
	src_2_dst_[ source_id ] = dst_2_src_.size();
}

platform::Size subset_mapping::source_size() const { return src_2_dst_.size(); }

platform::Size subset_mapping::destination_size() const { return dst_2_src_.size(); }

platform::Size subset_mapping::s2d( platform::Size source_id ) const { return src_2_dst_[ source_id ]; }

platform::Size subset_mapping::d2s( platform::Size destination_id ) const { return dst_2_src_[ destination_id ]; }

bool subset_mapping::source_id_is_mapped( platform::Size source_id ) const { return src_2_dst_[ source_id ] != UNMAPPED; }

}
