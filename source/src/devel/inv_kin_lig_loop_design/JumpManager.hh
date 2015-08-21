// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/InvKinLigLoopDesign/JumpManager.hh
///
/// @brief
/// @author

#ifndef DEVEL_INVKINLIGLOOPDESIGN_JUMPMANAGER_HH
#define DEVEL_INVKINLIGLOOPDESIGN_JUMPMANAGER_HH


#include <core/pose/Pose.hh>
#include <core/kinematics/Jump.hh>
#include <devel/inv_kin_lig_loop_design/Loop.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace devel {

namespace inv_kin_lig_loop_design {

using platform::Size;

struct JumpManager {

	// this function is required because it is not allowed because
	// not all atoms are allowed to be jump atoms. it returns the
	// closest atom which is allowed to be in a jump.
	bool is_allowable_jump_atom( core::pose::Pose& pose, int seqpos, int atom );
	string get_jump_atom_for_hbond( core::pose::Pose& pose, int seqpos, string hbond_atom );

	// this function sets a jump that when applied to the jump
	// between seqpos_from.jump_atom_from and seqpos_to.jump_atom_to
	// causes there to be a hydrogen bond between hbond_atom_from
	// and hbond_atom_to. NB: hbond_atom_from and to should be the
	// acceptor and hydrogen atoms irrespectively
	core::kinematics::RT get_hbond_rt( TagCOP tag, bool h2a );
	void set_random_hbond_jump(core::pose::Pose& pose, Loop const& loop );

	//core::kinematics::Jump get_jump_from_template( TagCOP tag_segment );
	core::kinematics::RT get_template_rt( TagCOP tag_segment );
	void set_template_jump( core::pose::Pose& pose, Loop const& loop );

}; // struct JumpManager

} // namespace LoopDesign

} // namespace Devel

#endif // DEVEL_LOOPDESIGN_JUMPMANAGER_HH
