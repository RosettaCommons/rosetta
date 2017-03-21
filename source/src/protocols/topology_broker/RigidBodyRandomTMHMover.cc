// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/topology_broker/RigidBodyRandomTMHMover.hh>
#include <protocols/topology_broker/TopologyClaimer.hh>

// Package headers
#include <protocols/rigid/RB_geometry.hh>
#include <core/pose/PDBInfo.hh>
// Rosetta Headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility>
#include <utility/tag/Tag.hh>

// Random number generator
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/TorsionID_Range.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/xyz.functions.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <string>
#include <core/id/TorsionID_Range.hh>
#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/tag/Tag.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace topology_broker {


static THREAD_LOCAL basic::Tracer TR( "protocols.moves.RigidBodyMover" );
static THREAD_LOCAL basic::Tracer TRBM( "protocols.moves.RigidBodyMover" );

RigidBodyRandomTMHMover::RigidBodyRandomTMHMover(){}
RigidBodyRandomTMHMover::~RigidBodyRandomTMHMover()= default;
RigidBodyRandomTMHMover::RigidBodyRandomTMHMover(core::Real max_trans, core::Real rotation_mag, core::Real translation_mag, core::Size tmhelix,
	protocols::topology_broker::TopologyClaimerOP claimer)
:    max_trans_(max_trans),

	trans_mag_in_(translation_mag),
	rot_mag_in_(rotation_mag),
	num_jump_(tmhelix),
	claimer_(std::move(claimer))
{
	if ( TR.Trace.visible() ) {
		TR.Trace << "max_trans:  " << max_trans_ << std::endl;
		TR.Trace << "rot_mag_in:  " << rot_mag_in_ << std::endl;
		TR.Trace << "trans_mag_in:  " << trans_mag_in_ << std::endl;
		TR.Trace << "num_jump:  " << num_jump_ << std::endl;
	}
}

void
RigidBodyRandomTMHMover::apply(core::pose::Pose& pose)
{
	//Choose a random jump
	core::Size random_jump_num = static_cast<core::Size>(numeric::random::rg().random_range(1,num_jump_));
	if ( TR.Trace.visible() ) { TR.Trace << "random_jump chosen:  " << random_jump_num << " " << pose.fold_tree().jump_edge(random_jump_num).start() << " " <<
		pose.fold_tree().jump_edge(random_jump_num).stop() << std::endl;
	}

	//Get the vector for this helix CoM (it is the second residue of the jump, the first is the virt res)
	core::Size current_span_CoM(pose.fold_tree().jump_edge(random_jump_num).stop());
	core::Vector current_rb_centroid(pose.residue(current_span_CoM).xyz("CA"));
	core::Vector original_centroid(0.0,0.0,0.0);
	core::pose::PoseOP claimer_pose;

	//If the claimer does have a PoseOP, get it
	if ( TR.Trace.visible() ) {
		TR.Trace << "Claimer is:  " << claimer_->type() << std::endl;
	}
	//get_pose_from_claimer() returns a PoseOP; NULL pointer in base class
	if ( claimer_->get_pose_from_claimer() ) {
		claimer_pose = claimer_->get_pose_from_claimer();
	} else {
		utility_exit_with_message("At least one TopologyClaimer needs to have get_pose_from_claimer() implemented for this mover to work");
	}

	//assert if claimer_pose doesn't point to anything at this point
	debug_assert(claimer_pose);

	//Get the claimer_pose's CoM for this helix.  This is the position to which we want to move to if we need to
	core::Size claimer_pose_current_span_CoM = claimer_pose->fold_tree().jump_edge(random_jump_num).stop();

	//The current pose and the claimer pose CoM for this helix should be the same!
	debug_assert(claimer_pose_current_span_CoM == current_span_CoM);

	//How far away are we now from the original starting position?
	original_centroid = claimer_pose->residue(claimer_pose_current_span_CoM).xyz("CA");
	if ( TR.Trace.visible() ) {
		TR.Trace << "original_centroid_xyz:  " << original_centroid.x() << " " << original_centroid.y() << " " << original_centroid.z() << std::endl;
	}

	core::Vector trans_vec = original_centroid - current_rb_centroid;
	core::Real trans_len = trans_vec.length();
	if ( TR.Trace.visible() ) {
		TR.Trace << "distance between current helix CoM and original position:  " << trans_len << std::endl;
	}

	//if we're too far away, just move back to original position from this claimer pose
	if ( trans_len>=max_trans_ ) {
		if ( TR.Trace.visible() ) {
			TR.Trace << "trans_len > max_trans.  Calling RigidBodyTransMover to move back to original position" << std::endl;
		}
		protocols::rigid::RigidBodyTransMover mover(pose, random_jump_num);
		mover.step_size(trans_len);
		mover.trans_axis(trans_vec);
		mover.apply(pose);
		current_rb_centroid = pose.residue(current_span_CoM).xyz("CA");
		if ( TR.Trace.visible() ) {
			TR.Trace << "current_rb_centroid position is now:  " << current_rb_centroid.x() << " " << current_rb_centroid.y() << " " << current_rb_centroid.z() << std::endl;
		}
		//if not too far away, just do RB perturbations as usual
	} else {
		if ( TR.Trace.visible() ) {
			TR.Trace << "trans_len <= max_trans.  Calling RigidBodyPerturbMover as usual" << std::endl;
		}
		rigid::RigidBodyPerturbMover mover(random_jump_num,rot_mag_in_,trans_mag_in_);
		mover.apply(pose);
	}
}

std::string
RigidBodyRandomTMHMover::get_name() const
{
	return "RigidBodyRandomTMHMover";
}

}
}

