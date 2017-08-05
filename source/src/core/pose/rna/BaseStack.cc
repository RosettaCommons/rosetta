// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/BaseStack.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/rna/BaseStack.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "core.pose.rna.BaseStack" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace rna {

using namespace core::chemical::rna;

///////////////////////////////////////////////////////////////////
BaseStack::BaseStack():
	res1_( 0 ),
	res2_( 0 ),
	orientation_( ANY_BASE_DOUBLET_ORIENTATION ),
	which_side_( ANY_BASE_STACK_SIDE )
{
}

BaseStack
BaseStack::flipped() const
{
	BaseStackWhichSide flipped_which_side( ANY_BASE_STACK_SIDE );
	if ( which_side_ == ABOVE ) {
		if ( orientation_ == ANTIPARALLEL ) {
			flipped_which_side = ABOVE;
		} else if ( orientation_ == PARALLEL ) {
			flipped_which_side = BELOW;
		}
	} else if ( which_side_ == BELOW ) {
		if ( orientation_ == ANTIPARALLEL ) {
			flipped_which_side = BELOW;
		} else if ( orientation_ == PARALLEL ) {
			flipped_which_side = ABOVE;
		}
	}
	return BaseStack( res2_, res1_, orientation_,
		flipped_which_side );
}

///////////////////////////////////////////////////////////////////
bool
operator < ( BaseStack const & lhs, BaseStack const & rhs ){
	// Priorities are res1_, then res2_, then orientation_, then which_side_
	if ( lhs.res1_ < rhs.res1_ ) {
		return true;
	}

	if ( lhs.res1_ == rhs.res1_ ) {
		if ( lhs.res2_ < rhs.res2_ ) {
			return true;
		}

		if ( lhs.res2_ == rhs.res2_ ) {
			if ( lhs.orientation_ < rhs.orientation_ ) {
				return true;
			}
			if ( lhs.orientation_ == rhs.orientation_ ) {
				return ( lhs.which_side_ < rhs.which_side_ );
			}
		}
	}
	return false;
}


///////////////////////////////////////////////////////////////////
std::ostream &
operator << ( std::ostream & out, BaseStack const & s )
{
	out << s.res1_ << " " << s.res2_ << " "
		<< get_full_orientation_from_num( s.orientation_ ) << " " << get_full_side_from_num( s.which_side_ );
	return out;
}

} //rna
} //pose
} //core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::rna::BaseStack::save( Archive & arc ) const {
	arc( CEREAL_NVP( res1_ ) ); // Size
	arc( CEREAL_NVP( res2_ ) ); // Size
	arc( CEREAL_NVP( orientation_ ) ); // enum core::chemical::rna::BaseDoubletOrientation
	arc( CEREAL_NVP( which_side_ ) ); // enum core::chemical::rna::BaseStackWhichSide
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::rna::BaseStack::load( Archive & arc ) {
	arc( res1_ ); // Size
	arc( res2_ ); // Size
	arc( orientation_ ); // enum core::chemical::rna::BaseDoubletOrientation
	arc( which_side_ ); // enum core::chemical::rna::BaseStackWhichSide
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::rna::BaseStack );
CEREAL_REGISTER_TYPE( core::pose::rna::BaseStack )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_rna_BaseStack )
#endif // SERIALIZATION
