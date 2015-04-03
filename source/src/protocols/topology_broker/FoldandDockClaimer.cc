// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FoldandDockClaimer
/// @brief Fold-and-dock
/// @author Ingemar Andre

// Unit Headers
#include <protocols/topology_broker/FoldandDockClaimer.hh>
#include <protocols/symmetric_docking/SymFoldandDockRbTrialMover.hh>
#include <protocols/symmetric_docking/SymFoldandDockSlideTrialMover.hh>
#include <protocols/simple_moves/symmetry/SymDockingInitialPerturbation.hh>
#include <protocols/symmetric_docking/SymFoldandDockMoveRbJumpMover.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/CutClaim.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <basic/Tracer.hh>

// Utility header
#include <core/pose/symmetry/util.hh>

#include <numeric/random/random.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/fold_and_dock.OptionKeys.gen.hh>


#include <core/chemical/AtomType.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/moves/MoverContainer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>

//C++ Headers
#include <utility>

// Project Headers

static thread_local basic::Tracer tr( "protocols.topo_broker.fold_and_dock", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace core;


FoldandDockClaimer::FoldandDockClaimer() {}

FoldandDockClaimer::FoldandDockClaimer( pose::Pose const& input_pose ) :
	input_pose_(input_pose)
{}

	//clone
TopologyClaimerOP
FoldandDockClaimer::clone() const {
	return TopologyClaimerOP( new FoldandDockClaimer( *this ) );
}

/// @brief type() is specifying the output name of the TopologyClaimer
std::string
FoldandDockClaimer::type() const {
	return _static_type_name();
}

std::string
FoldandDockClaimer::_static_type_name() {
	return "FoldandDockClaimer";
}

void
FoldandDockClaimer::add_mover(
    moves::RandomMover& random_mover,
		core::pose::Pose const& /*pose*/,
		abinitio::StageID stageID,  /*abinitio sampler stage */
		core::scoring::ScoreFunction const& scorefxn,
		core::Real /*progress  progress within stage */
)
{
	using namespace basic::options;

	moves::MoverOP move_anchor_mover( new symmetric_docking::SymFoldandDockMoveRbJumpMover );
	moves::MoverOP rb_trial_mover(
		(stageID==abinitio::STAGE_4) ?
		new symmetric_docking::SymFoldandDockRbTrialMover( scorefxn.get_self_ptr(), true ) :
		new symmetric_docking::SymFoldandDockRbTrialMover( scorefxn.get_self_ptr() )  // smooth RB moves in stage 4
	);
	moves::MoverOP slide_mover( new symmetric_docking::SymFoldandDockSlideTrialMover );
	core::Real move_anchor_weight(option[ OptionKeys::fold_and_dock::move_anchor_frequency ]()),
	           rb_weight(option[ OptionKeys::fold_and_dock::rigid_body_frequency ]()),
	           slide_weight(option[ OptionKeys::fold_and_dock::slide_contact_frequency ]());

	if (move_anchor_weight > 0) random_mover.add_mover( move_anchor_mover, move_anchor_weight );
	random_mover.add_mover( rb_trial_mover, rb_weight );
	random_mover.add_mover( slide_mover, slide_weight );
}

void FoldandDockClaimer::initialize_dofs(
	core::pose::Pose& pose,
	claims::DofClaims const& init_dofs,
	claims::DofClaims& /*failed_to_init*/ ) {

	using namespace core::conformation::symmetry;

	// Setup symmetry if we have nit already done it
	// slide chains into contact
	protocols::simple_moves::symmetry::SetupForSymmetryMoverOP setup_mover( new
		protocols::simple_moves::symmetry::SetupForSymmetryMover );
	setup_mover->slide_into_contact(true);
	if ( !core::pose::symmetry::is_symmetric( pose ) ) {
		setup_mover->apply( pose ); // calls SymDockingInitialPerturbation
		assert( core::pose::symmetry::is_symmetric( pose ) );
		// Save the pose into input pose
		input_pose_ = pose;
	} else {
		input_pose_ = pose;
		// Randomize the rigid body
		protocols::simple_moves::symmetry::SymDockingInitialPerturbation initial( true /*slide into contact*/ );
		initial.apply( pose );
	}

	// Setup the movemap
	//SymmetricConformation & symm_conf (dynamic_cast<SymmetricConformation & > ( pose.conformation()) );
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

void FoldandDockClaimer::generate_claims( claims::DofClaims& new_claims ) {
	// Set all cuts to real cuts. We don't want to close any of them...
	utility::vector1< int > cuts( input_pose_.conformation().fold_tree().cutpoints() );
	for ( Size i = 1; i <= cuts.size(); ++i ) {
		new_claims.push_back( claims::DofClaimOP( new claims::CutClaim( get_self_weak_ptr(), std::make_pair( Parent::label(), cuts[i]),
																								claims::DofClaim::INIT /* for now... eventually CAN_INIT ? */ ) ) );
	}
}


} //topology_broker
} //protocols
