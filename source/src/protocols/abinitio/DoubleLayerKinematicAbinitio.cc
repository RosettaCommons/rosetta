// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file KinematicTaskCenter
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

// Unit Headers
#include <protocols/abinitio/DoubleLayerKinematicAbinitio.hh>

// Package Headers
#include <protocols/abinitio/ResolutionSwitcher.hh>
#include <protocols/simple_moves/FragmentMover.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <core/kinematics/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/loops/Loop.hh>


#include <basic/options/option.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <numeric/random/random.hh>
#include <basic/Tracer.hh>


//// C++ headers


// option key includes

#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.general_abinitio", basic::t_info );

namespace protocols {
namespace abinitio {

using namespace core;

DoubleLayerKinematicAbinitio::~DoubleLayerKinematicAbinitio() = default;

std::string
DoubleLayerKinematicAbinitio::get_name() const {
	return "DoubleLayerKinematicAbinitio";
}

///////////////////////////////////////////////////////////////////////
/// @brief Select loop set at random using skip rate
void DoubleLayerKinematicAbinitio::select_core_loops(
	loops::Loops& loops_out
) const
{
	loops_out.clear();
	int ntries = 0;
	while ( loops_out.size() == 0 && ntries++ < 50 ) {
		 for ( auto const & rigid_loop : rigid_loops_ ) {
			if ( numeric::random::rg().uniform() >= rigid_loop.skip_rate() )  {
				loops_out.push_back( rigid_loop );
			}
		}
	}
	if ( loops_out.size() == 0 ) {
		loops_out = rigid_loops_;
	}
} // void LoopRebuild::select_loops


//@brief create a new fold-tree and movemap --- a KinematicControl object
// a basic generalized protocol:
// loops are determined: if loops present only loops move, otherwise everything moves
// get jumps from JumpDef ( as in JumpingFoldCst, for instance beta-sheet jumps )
// combine jumps with additional jumps that keep the rigid portions together ( looprlx-type )
// set movemap to allow only sampling of flexible parts. call sampling protocol
KinematicControlOP DoubleLayerKinematicAbinitio::new_kinematics( pose::Pose &pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// select loops
	loops::Loops rigid_core;
	if ( rigid_loops_.size() > 0 ) {
		select_core_loops( rigid_core );
	}
	tr.Debug << rigid_core << std::endl;

	KinematicControlOP current_kinematics( nullptr );
	if ( rigid_core.size() && coordinate_constraint_weight_ > 0.0 ) {
		current_kinematics = KinematicControlOP( new CoordinateConstraintKC( false /*ramp*/, coordinate_constraint_weight_ ) );
	} else {
		current_kinematics = KinematicControlOP( new KinematicControl );
	}
	loops::Loops flexible_part( rigid_core.invert( pose.total_residue() ) );

	bool loop_file_is_present = option[ OptionKeys::loops::mm_loop_file ].user();
	loops::Loops mmloops( loop_file_is_present ) ;

	if ( !loop_file_is_present ) {
		mmloops = flexible_part;
	}

	//figure out movemap!
	kinematics::MoveMapOP movemap( new kinematics::MoveMap );
	movemap->set_jump( true ); //why is that here ?

	if ( mmloops.size() && coordinate_constraint_weight_ == 0.0 ) {
		rigid_core.switch_movemap( *movemap, id::BB, false );
		mmloops.switch_movemap( *movemap, id::BB, true );
		flexible_part.switch_movemap( *movemap, id::BB, true );
	} else {
		movemap->set_bb( true );
	}


	//figure out flexible jumps - jump-movers
	bool success( true );
	current_kinematics->set_final_fold_tree( pose.fold_tree() );
	success = parse_jump_def( current_kinematics, movemap );

	current_kinematics->set_movemap( movemap );
	// change fold-tree such that only loop parts can move
	if ( rigid_core.size() ) {
		success &= add_rigidity_jumps( rigid_core, current_kinematics );
		if ( !success ) {
			tr.Warning << "[WARNING] was not able to fix rigid regions with jumps...retry" << std::endl;
			return nullptr;
		}
	}

	pose.fold_tree( current_kinematics->sampling_fold_tree() );

	set_extended_torsions_and_idealize_loops( pose, extended_loops_ );

	if ( rigid_core.size() &&  coordinate_constraint_weight_ != 0.0 ) {
		/*success = */add_coord_cst( rigid_core, pose );
	}
	return current_kinematics;
}


//@brief sampling: simple version: get new kinematics (movemap+jumps) and call sampling protocol.
//overwrite this guy if you want to do more stuff ... i.e., extend loops if things didn't work out in the first place.
bool  DoubleLayerKinematicAbinitio::inner_loop( core::pose::Pose& pose ) {
	bool success( false );

	Size fail( 0 );
	current_kinematics_ = nullptr;
	while ( fail++ <= 10 && !current_kinematics_ ) {// get new setup
		//this may add constraints to the pose ...or should this be handled via the KinematicControl object?!
		current_kinematics_ = new_kinematics( pose );
	}

	//debug output
	if ( current_kinematics_ && tr.Info.visible() ) {
		tr.Info << "kinematic choice:\n";
		core::kinematics::simple_visualize_fold_tree_and_movemap(
			current_kinematics_->sampling_fold_tree(),
			current_kinematics_->movemap(),
			tr.Info );
		tr.Info << "\n" << jumping::JumpSample( current_kinematics_->sampling_fold_tree() );
		tr.Info << "\nfinal_fold-tree:\n";
		core::kinematics::simple_visualize_fold_tree( current_kinematics_->final_fold_tree(), tr.Info );
	}

	// if setup valid...
	if ( current_kinematics_ ) {
		// sample with this setup

		// first sample only extended in stage1
		kinematics::MoveMapOP extended_movemap( new kinematics::MoveMap( current_kinematics()->movemap() ) );
		extended_movemap->set_bb( false );
		extended_loops_.switch_movemap( *extended_movemap, id::BB, true );

		KinematicControlOP stage1_kinematics( new KinematicControl( *current_kinematics() ) );
		stage1_kinematics->set_movemap( extended_movemap );
		stage1_sampler_->set_kinematics( stage1_kinematics );
		stage1_sampler_->set_current_tag( get_current_tag() + "_presampled" );
		stage1_sampler_->apply( pose );

		sampling_protocol()->set_kinematics( current_kinematics() );
		sampling_protocol()->apply( pose );
		success = ( sampling_protocol()->get_last_move_status() == moves::MS_SUCCESS );
	}
	// done...
	return success;
}

}
}
