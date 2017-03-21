// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/RelativeConnectRight.cc
/// @brief version of ConnectRight instruction that depends upon results from
///  another BuildInstruction
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/build/RelativeConnectRight.hh>
#include <protocols/forge/build/RelativeSequencePosition.hh>

// project headers
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>


#ifdef WIN32
// apparently this is required for a Visual Studio build.
#include <core/conformation/Residue.hh>
#endif


namespace protocols {
namespace forge {
namespace build {


/// @brief default constructor
RelativeConnectRight::RelativeConnectRight() :
	Super()
{}


/// @brief RelativeSequencePosition + position on-right jump constructor
/// @param[in] rp RelativeSequencePosition defining the type of computation to perform.
///  (will be cloned)
/// @param[in] right_position connect at this position on 'pose_right'
/// @param[in] pose_right connect this pose to the right of pose_left when
///  modify( pose_left ) is called
RelativeConnectRight::RelativeConnectRight(
	RelativeSequencePositionOP const & rp,
	Size const right_position,
	Pose const & pose_right
) :
	Super( 0u, right_position, pose_right ),
	rp_( rp->clone() )
{}


/// @brief copy constructor
RelativeConnectRight::RelativeConnectRight( RelativeConnectRight const & rval ) :
	Super( rval ),
	rp_( rval.rp_->clone() )
{}


/// @brief default destructor
RelativeConnectRight::~RelativeConnectRight() {}


/// @brief copy assignment
RelativeConnectRight & RelativeConnectRight::operator =( RelativeConnectRight const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );
		rp_ = rval.rp_->clone();
	}

	return *this;
}


/// @brief clone this object
BuildInstructionOP RelativeConnectRight::clone() const {
	return BuildInstructionOP( new RelativeConnectRight( *this ) );
}


/// @brief do the actual work of modifying the Pose
void RelativeConnectRight::modify_impl( Pose & pose_left ) {
	// modify 'left_position' to refer to proper jump takeoff wrt RelativeSequencePosition
	// function object
	debug_assert( n_dependencies() == 1 );
	left_position( (*rp_)( *dependencies().begin() ) );

	// do the actual work
	Super::modify_impl( pose_left );
}


/// @brief return set of any fixed positions necessary with respect to the original
///  interval and original Pose numbering
/// @return always empty set, no fixed positions
/// @remarks Used for ensuring build regions for instructions do not overlap and
///  so that jumps may be placed correctly.  There is currently no way to
///  represent the dependent fixed position, so we're forced to return an empty
///  set.
RelativeConnectRight::Positions RelativeConnectRight::original_fixed_positions() const {
	return Positions();
}


} // namespace build
} // namespace forge
} // namespace protocols
