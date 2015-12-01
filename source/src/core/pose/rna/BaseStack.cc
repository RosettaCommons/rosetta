// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/rna/BaseStack.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/rna/BaseStack.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "core.pose.rna.BaseStack" );

namespace core {
namespace pose {
namespace rna {

///////////////////////////////////////////////////////////////////
BaseStack::BaseStack():
	res1_( 0 ),
	res2_( 0 ),
	orientation_( ANY_BASE_DOUBLET_ORIENTATION ),
	which_side_( ANY_BASE_STACK_SIDE )
{
}

///////////////////////////////////////////////////////////////////
bool
operator < ( BaseStack const & lhs, BaseStack const & rhs ){
	//There must be a more elegant way to do this...
	if ( lhs.res1_ < rhs.res1_ ) {
		return true;
	}  else if ( lhs.res1_ == rhs.res1_ ) {
		if ( lhs.res2_ < rhs.res2_ ) {
			return true;
		} else if ( lhs.res2_ == rhs.res2_ ) {
			if ( lhs.orientation_ < rhs.orientation_ ) {
				return true;
			} else if ( lhs.orientation_ == rhs.orientation_ ) {
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
	out << s.res1_ << " " << s.res2_ << " " <<  s.orientation_ << " " << s.which_side_;
	return out;
}

} //rna
} //pose
} //core
