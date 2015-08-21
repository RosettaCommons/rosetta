// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FibrilModelingClaimer
/// @brief Fibril Modeling
/// @author Lin Jiang

// Unit Headers
#include <protocols/topology_broker/FibrilModelingClaimer.hh>
#include <protocols/symmetric_docking/SymFoldandDockRbTrialMover.hh>
#include <protocols/symmetric_docking/SymFoldandDockSlideTrialMover.hh>
#include <protocols/simple_moves/symmetry/SymDockingInitialPerturbation.hh>
#include <protocols/symmetric_docking/SymFoldandDockMoveRbJumpMover.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/BBClaim.hh>
#include <protocols/topology_broker/claims/CutClaim.hh>
#include <protocols/fibril/SetupForFibrilMover.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/fold_and_dock.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility header
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <numeric/random/random.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/LoopsFileIO.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/moves/MoverContainer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>

static thread_local basic::Tracer tr( "protocols.topo_broker.fibril_modeling", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace core;


FibrilModelingClaimer::FibrilModelingClaimer() { bAlign_ = false; sequence_shift_ = 0;}

FibrilModelingClaimer::FibrilModelingClaimer( pose::Pose const& input_pose, loops::Loops rigid, int shift ):
	input_pose_(input_pose)
{
	rigid_core_ = rigid;
	input_rigid_core_ = rigid_core_ ;
	sequence_shift_ = shift;
	input_rigid_core_.make_sequence_shift(sequence_shift_);
	bAlign_ = true;
}

FibrilModelingClaimer::FibrilModelingClaimer( pose::Pose const& input_pose, loops::Loops rigid, loops::Loops input_rigid ):
	input_pose_(input_pose)
{
	rigid_core_ = rigid;
	input_rigid_core_ = input_rigid ;
	sequence_shift_ = 0;
	bAlign_ = true;
}

bool FibrilModelingClaimer::read_tag( std::string tag, std::istream& is ) {

	if ( tag == "pdb" || tag == "PDB" || tag == "pdb:" || tag == "PDB_FILE" ) {
		std::string file;
		is >> file;
		core::import_pose::pose_from_pdb( input_pose_,
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ),
			file );
		bAlign_ = true;
		runtime_assert( input_pose_.is_fullatom() );
	} else if ( tag == "sequence_shift" ) {
		is >> sequence_shift_ ;
	} else if ( tag == "REGION" ) {
		loops::PoseNumberedLoopFileReader reader;
		reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
		loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file( is, type(), false /*no strict checking */ );
		rigid_core_ = loops::Loops( loops );
		input_rigid_core_ = rigid_core_;
		input_rigid_core_.make_sequence_shift(sequence_shift_);
	} else if ( tag == "INPUT_REGION" ) {
		input_rigid_core_.clear();
		loops::PoseNumberedLoopFileReader loop_file_reader;
		loop_file_reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
		loops::SerializedLoopList loops = loop_file_reader.read_pose_numbered_loops_file( is, type(), false /*no strict checking */ );
		input_rigid_core_ = loops::Loops( loops );
	} else return Parent::read_tag( tag, is );
	return true;
}

void
FibrilModelingClaimer::make_fibril( pose::Pose & pose )
{
	using namespace core::conformation::symmetry;

	// Setup symmetry if we have not already done it
	if ( !core::pose::symmetry::is_symmetric( pose ) ) {
		protocols::fibril::SetupForFibrilMoverOP setup_mover( new protocols::fibril::SetupForFibrilMover );
		if ( bAlign_ ) {
			setup_mover->align( pose, input_pose_ , rigid_core_, input_rigid_core_ );
		}
		setup_mover->apply( pose );
	}
	assert( core::pose::symmetry::is_symmetric( pose ) );
	// Save the symmetry info
	// input_pose_ = pose;
	symminfo_ = core::pose::symmetry::symmetry_info(pose)->clone();

}

void
FibrilModelingClaimer::add_mover(
	moves::RandomMover& random_mover,
	core::pose::Pose const& /*pose*/,
	abinitio::StageID /*stageID,  abinitio sampler stage */,
	core::scoring::ScoreFunction const& scorefxn,
	core::Real /*progress  progress within stage */
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::fold_and_dock;

	core::Real move_anchor_weight(1.0), rb_weight, slide_weight;

	if ( option( move_anchor_points ).user() ) {
		moves::MoverOP move_anchor_mover( new symmetric_docking::SymFoldandDockMoveRbJumpMover );
		random_mover.add_mover( move_anchor_mover, move_anchor_weight );
	}

	rb_weight = option( rigid_body_frequency );
	moves::MoverOP rb_trial_mover( new symmetric_docking::SymFoldandDockRbTrialMover( scorefxn.get_self_ptr() ) );
	random_mover.add_mover( rb_trial_mover, rb_weight );

	slide_weight = option( slide_contact_frequency );
	moves::MoverOP slide_mover( new symmetric_docking::SymFoldandDockSlideTrialMover );
	random_mover.add_mover( slide_mover, slide_weight );

}

void FibrilModelingClaimer::initialize_dofs(
	core::pose::Pose& pose,
	claims::DofClaims const& init_dofs,
	claims::DofClaims& /*failed_to_init*/ ) {

	using namespace core::conformation::symmetry;

	// Setup symmetry if we have not already done it
	make_fibril( pose );

	// Randomize the rigid body
	protocols::simple_moves::symmetry::SymDockingInitialPerturbation initial( false /*slide into contact*/ );
	initial.apply( pose );
	// Setup the movemap
	//SymmetricConformation & symm_conf (
	//     dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

	kinematics::MoveMapOP movemap( new kinematics::MoveMap() );
	movemap->set_bb( true );
	movemap->set_jump( false );
	core::pose::symmetry::make_symmetric_movemap( pose, *movemap );

	for ( claims::DofClaims::const_iterator it = init_dofs.begin(), eit = init_dofs.end();
			it != eit; ++it ) {
		if ( (*it)->owner().lock().get() == this ) {
			(*it)->toggle( *movemap, true );
		}
	}
}

void FibrilModelingClaimer::generate_claims( claims::DofClaims& new_claims ) {
	// Set all cuts to real cuts. We don't want to close any of them...
	utility::vector1< int > cuts( input_pose_.conformation().fold_tree().cutpoints() );
	for ( Size i = 1; i <= cuts.size(); ++i ) {
		new_claims.push_back( claims::DofClaimOP( new claims::CutClaim( get_self_weak_ptr(), std::make_pair( Parent::label(), cuts[i]),
			claims::DofClaim::INIT /* for now... eventually CAN_INIT ? */ ) ) );
	}
}


bool FibrilModelingClaimer::allow_claim( claims::DofClaim const& foreign_claim ) {

	using namespace core::conformation::symmetry;

	if ( foreign_claim.owner().lock().get() == this ) return true; // always allow your own claims!

	//std::cout<<"claim_pos "<<foreign_claim.pos( 1 )<<" "<<is_symmetric( input_pose_ )<<std::endl;
	if ( core::conformation::symmetry::is_symmetric( *symminfo_ ) ) {


		// check foreign claim
		claims::BBClaim const *bb_ptr( dynamic_cast< const claims::BBClaim* >( &foreign_claim ) );
		if ( bb_ptr ) {
			if ( ! symminfo_->bb_is_independent( bb_ptr->global_position() ) ) return false;
		} // DofClaim::BB

		//if ( foreign_claim.type() == DofClaim::JUMP ) {
		// return false;
		//}

		//if ( foreign_claim.type() == DofClaim::CUT) {
		// if( ! input_symminfo_.bb_is_independent( foreign_claim.pos( 1 ) ) ) return false;
		//}
	}

	return true;
} // FibrilModelingClaimer::allow_claim()


} //topology_broker
} //protocols
