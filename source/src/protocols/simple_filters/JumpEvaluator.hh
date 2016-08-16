// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file JumpEvaluator.hh
/// @brief
/// @details
///
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_simple_filters_JumpEvaluator_hh
#define INCLUDED_protocols_simple_filters_JumpEvaluator_hh


// Unit Headers

// Package Headers
#include <protocols/evaluation/PoseEvaluator.hh>

// Project Headers
#include <core/id/AtomID.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/Stub.hh>

#include <core/io/silent/silent.fwd.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


//// C++ headers

namespace protocols {
namespace simple_filters {

class JumpEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	JumpEvaluator( core::pose::Pose const& native_pose, Size jump_nr );

	virtual core::Real apply( core::pose::Pose& pose  ) const;

private:
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Size jump_nr_;
	core::id::AtomID up_jump_atom_;
	core::id::AtomID down_jump_atom_;
	core::id::StubID down_stub_;
	core::id::StubID up_stub_;
	core::kinematics::Stub native_up_;
	core::kinematics::Stub native_down_;
};

//@brief yields a column with the number of jumps in the pose
class JumpNrEvaluator : public evaluation::SingleValuePoseEvaluator< core::Size > {
public:
	JumpNrEvaluator() : evaluation::SingleValuePoseEvaluator< core::Size >( "nrjumps" ) {};
	virtual core::Size apply( core::pose::Pose& pose  ) const;
private:
};


}
}

#endif
