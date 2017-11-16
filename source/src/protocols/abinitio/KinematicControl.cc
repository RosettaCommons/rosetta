// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file KinematicControl.cc
/// @brief
/// @details
/// @author Oliver Lange


// Unit Headers
#include <protocols/abinitio/KinematicControl.hh>

// Package Headers

// Project Headers
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>


#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>


#include <core/scoring/ScoreFunction.hh>

#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/jumping/util.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility headers

#include <basic/Tracer.hh>


//// C++ headers

//Auto Headers
#include <core/chemical/AtomType.hh>
#include <core/pose/variant_util.hh>
#include <utility/vector1.hh>


static basic::Tracer tr( "protocols.abinitio", basic::t_info );

namespace protocols {
namespace abinitio {

using namespace core;

bool KinematicControl::prepare_pose_for_sampling( pose::Pose& pose ) const {
	pose.fold_tree( sampling_fold_tree() );
	if ( jump_mover() ) jump_mover()->apply_at_all_positions( pose ); //make sure each jump is initialized
	jumping::safe_secstruct( pose ); //make sure that secstruct is valid (in the sense of FragmentMover::valid_ss)
	return true;
}

//@brief find all cutpoints that are only present in the "sampling" fold-tree.
void
find_sampling_cuts(
	kinematics::FoldTree const& sampling,
	kinematics::FoldTree const& final,
	utility::vector1< Size >& sample_cuts )
{
	sample_cuts.clear();
	for ( Size i = 1; i <= (Size) sampling.num_cutpoint(); i++ ) {
		if ( !final.is_cutpoint( sampling.cutpoint( i ) ) ) sample_cuts.push_back( sampling.cutpoint( i ));
	}
}

void
KinematicControl::add_chainbreak_variants( pose::Pose &pose ) const {
	utility::vector1< Size > sample_cuts;
	find_sampling_cuts( pose.fold_tree(), final_fold_tree(), sample_cuts );
	for ( utility::vector1< Size >::const_iterator it = sample_cuts.begin(), eit = sample_cuts.end();
			it != eit; ++ it ) {
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, *it );
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, *it+1 );
	}
}

void
KinematicControl::add_chainbreak_variants( pose::Pose &pose, Size max_dist, core::kinematics::ShortestPathInFoldTree const& sp) const {
	//remove_chainbreaks( pose ); not necessary if max_dist is monotonoically increaseing
	utility::vector1< Size > sample_cuts;
	find_sampling_cuts( pose.fold_tree(), final_fold_tree(), sample_cuts );
	for ( utility::vector1< Size >::const_iterator it = sample_cuts.begin(), eit = sample_cuts.end();
			it != eit; ++ it ) {
		if ( sp.dist( *it, *it+1 ) <= max_dist ) {
			tr.Debug << "add chainbreak variant to residues " << *it << " and " << *it+1 << std::endl;
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, *it );
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, *it+1 );
		}
	}
}

KinematicControl::KinematicControl() {}
KinematicControl::~KinematicControl() = default;

void KinematicControl::set_movemap( core::kinematics::MoveMapCOP mm ) {
	movemap_ = mm;
	if ( jump_mover_ ) jump_mover_->set_movemap( movemap_ptr() );
}

void KinematicControl::set_strict_movemap( core::kinematics::MoveMapCOP mm ) {
	strict_movemap_ = mm;
}

core::kinematics::MoveMapCOP
KinematicControl::movemap_ptr() const {
	return movemap_;
}

core::kinematics::MoveMap const &
KinematicControl::movemap() const {
	return *movemap_;
}

//return a jump-Mover for jumps that you want to be sampled
simple_moves::FragmentMoverOP
KinematicControl::jump_mover() const {
	return jump_mover_;
}

//return a jump-Mover for jumps that you want to be sampled
void
KinematicControl::set_jump_mover( simple_moves::FragmentMoverOP jm ) {
	jump_mover_ = jm;
	if ( jump_mover_ && movemap_ ) jump_mover_->set_movemap( movemap_ );
}


void CoordinateConstraintKC::add_score_weights( scoring::ScoreFunction& scorefxn, Real progress ) const {
	scorefxn.set_weight( scoring::coordinate_constraint, ramp_ ? progress*final_weight_ : final_weight_ );
}

}
}
