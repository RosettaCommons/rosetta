// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/rigid/RotateJumpAxisMover.hh
/// @brief RotateJumpAxisMover header
/// @author Steven Lewis

#ifndef INCLUDED_protocols_rigid_RotateJumpAxisMover_hh
#define INCLUDED_protocols_rigid_RotateJumpAxisMover_hh

// Unit Headers
#include <protocols/rigid/RotateJumpAxisMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace rigid {

/// @details This mover rotates a jump transform.  Its original use was to rotate a freely rotateable zinc-histidine bond as emulated by an atom-to-atom fixed-length jump.  It recalculates the Stubs for the jump and applies the new jump, resulting in an N degree rotation of one partner about the axis between the histidine nitrogen and the zinc.  It will work for any jump but is intended for atom-atom jumps (not residue-residue jumps).  It will choose an angle from the uniform random distribution bounded by inputs (defaults to (-180,180]); if you want a particular value then set the limits equal.
class RotateJumpAxisMover : public protocols::moves::Mover {

public:

	/// @brief constructor for random distribution (just needs rb_jump_num)
	RotateJumpAxisMover( core::Size const rb_jump_num );

	/// @brief constructor for range - these angles are in degrees, not radians!
	RotateJumpAxisMover( core::Size const rb_jump_num, core::Angle const upper, core::Angle const lower );

	/// @brief constructor for single value - these angles are in degrees, not radians!
	RotateJumpAxisMover( core::Size const rb_jump_num, core::Angle const angle );

	~RotateJumpAxisMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:
	core::Angle calc_angle();

	core::Size const rb_jump_num_;
	core::Angle const upper_angle_; //these angles are in degrees, not radians!
	core::Angle const lower_angle_; //these angles are in degrees, not radians!

};//end RotateJumpAxisMover

}//namespace rigid
}//namespace protocols

#endif // INCLUDED_protocols_rigid_RotateJumpAxisMover_HH
