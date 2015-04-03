// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file KinematicControl
/// @brief  this class will be handled to a SampleProtocol as a control instance
/// @details responsibilities:
///           know which chainbreaks to penalize and close
///           know which jumps to use during sampling, which (if any) to keep after loop-closing
///           supply a JumpMover if jumps should be moved
///           supply a MoveMap
///           supply a "StrictMoveMap": the protocol should not move anything that is dissallowed in strict_movemap(),
///                      it should try to move just stuff in movemap()
/// should this class also know how to ramp score terms ?
/// handle the titration of constraints ?
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_KinematicControl_hh
#define INCLUDED_protocols_abinitio_KinematicControl_hh

// Unit Headers
#include <protocols/abinitio/KinematicControl.fwd.hh>

// Package Headers

// Project Headers
#include <core/types.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


#include <protocols/simple_moves/FragmentMover.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

// enum JumpType {
//   BETA = 1,
//   HELIX,
//   LOOPRLX //never sample orientation on these guys
// };

// class AnnotatedJump {
// public:
//   core::Size start_;
//   core::Size end_;
//   bool bPermanent_;
//   JumpType type_; // does this makes sense ?
// };


class KinematicControl : public utility::pointer::ReferenceCount {
public:
	KinematicControl();
	virtual ~KinematicControl();

	//@brief setup things in pose: e.g., set correct fold-tree, set jump-geometries to initial values...
	bool prepare_pose_for_sampling( core::pose::Pose& pose ) const;

	void add_chainbreak_variants( core::pose::Pose& pose )  const;

	//@brief switch on chainbreaks-scores only if the separation in the foldtree is smaller than the threshold <max_dist>
	void
	add_chainbreak_variants( core::pose::Pose &pose, core::Size max_dist, core::kinematics::ShortestPathInFoldTree const& ) const;

	//void remove_chainbreak_variants( core::pose::Pose& pose ) const;

	core::kinematics::FoldTree const& sampling_fold_tree() const {
		return sampling_fold_tree_;
	}

	core::kinematics::FoldTree const& final_fold_tree() const {
		return final_fold_tree_;
	}

	void set_sampling_fold_tree( core::kinematics::FoldTree const& f) {
		sampling_fold_tree_ = f;
	}
  //  void set_sample_cuts();//???
  //void set_sample_jumps(); //???

	void set_final_fold_tree( core::kinematics::FoldTree const& f) {
		final_fold_tree_ = f;
	}
	//void set_final_cuts();//???
	//void set_final_jumps(); //???

	void set_movemap( core::kinematics::MoveMapCOP mm );

	void set_strict_movemap( core::kinematics::MoveMapCOP mm );

	core::kinematics::MoveMapCOP movemap_ptr() const;

	core::kinematics::MoveMap const& movemap() const;

	//return a jump-Mover for jumps that you want to be sampled
	simple_moves::FragmentMoverOP jump_mover() const;

	//return a jump-Mover for jumps that you want to be sampled
	void set_jump_mover( simple_moves::FragmentMoverOP jm );

	virtual void add_score_weights( core::scoring::ScoreFunction&, core::Real /*progress*/ ) const {};

// 	loops::Loops const& rigid_zones() const {
// 		return rigid_;
// 	}

// 	void rigid_zones( loops::Loops const& rigid ) {
// 		rigid_ = rigid;
// 	}

private:
  core::kinematics::MoveMapCOP strict_movemap_;
  core::kinematics::MoveMapCOP movemap_;

  core::kinematics::FoldTree sampling_fold_tree_;
  core::kinematics::FoldTree final_fold_tree_;

  simple_moves::FragmentMoverOP jump_mover_;

// 	loops::Loops rigid_;
};

class CoordinateConstraintKC : public KinematicControl {
public:
	CoordinateConstraintKC( bool ramp, core::Real final_weight ) :
		ramp_ ( ramp ),
		final_weight_( final_weight )
	{}

	virtual void add_score_weights( core::scoring::ScoreFunction&, core::Real progress ) const;

private:
	bool ramp_;
	core::Real final_weight_;
};

}
}
#endif
