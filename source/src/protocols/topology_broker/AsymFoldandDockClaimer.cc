// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AsymFoldandDockClaimer
/// @brief Fold-and-dock
/// @author Ingemar Andre

// Unit Headers
#include <protocols/topology_broker/AsymFoldandDockClaimer.hh>
#include <protocols/simple_moves/asym_fold_and_dock/AsymFoldandDockRbTrialMover.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/simple_moves/asym_fold_and_dock/AsymFoldandDockMoveRbJumpMover.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/CutClaim.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/docking/util.hh>

// Utility headers
#include <numeric/random/random.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/fold_and_dock.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/moves/MoverContainer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//C++ Headers
#include <utility>

// Project Headers

static THREAD_LOCAL basic::Tracer tr( "protocols.topo_broker.asym_fold_and_dock", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace core;


AsymFoldandDockClaimer::AsymFoldandDockClaimer() {}

AsymFoldandDockClaimer::AsymFoldandDockClaimer( pose::Pose const& input_pose ) :
	input_pose_(input_pose)
{
	docking_local_refine_ = false;
	chain_break_res_ = 0;
}

//clone
TopologyClaimerOP
AsymFoldandDockClaimer::clone() const {
	return TopologyClaimerOP( new AsymFoldandDockClaimer( *this ) );
}

/// @brief type() is specifying the output name of the TopologyClaimer
std::string
AsymFoldandDockClaimer::type() const {
	return _static_type_name();
}

std::string
AsymFoldandDockClaimer::_static_type_name() {
	return "AsymFoldandDockClaimer";
}

void
AsymFoldandDockClaimer::add_mover(
	moves::RandomMover& random_mover,
	core::pose::Pose const& /*pose*/,
	abinitio::StageID stageID,  /*abinitio sampler stage */
	core::scoring::ScoreFunction const& scorefxn,
	core::Real /*progress  progress within stage */
)
{
	using namespace basic::options;
	using namespace protocols::loops;
	using namespace protocols::simple_moves;

	moves::MoverOP move_anchor_mover( new simple_moves::asym_fold_and_dock::AsymFoldandDockMoveRbJumpMover( chain_break_res_ ) );
	moves::MoverOP rb_trial_mover (
		(stageID==abinitio::STAGE_4) ?
		new simple_moves::asym_fold_and_dock::AsymFoldandDockRbTrialMover( scorefxn.get_self_ptr(), true ) :
		new simple_moves::asym_fold_and_dock::AsymFoldandDockRbTrialMover( scorefxn.get_self_ptr() )  // smooth RB moves in stage 4
	);
	moves::MoverOP slide_mover( new protocols::docking::DockingSlideIntoContact( docking_jump_ ) );
	core::Real move_anchor_weight(1.0),
		rb_weight(option[ OptionKeys::fold_and_dock::rigid_body_frequency ]()),
		slide_weight(option[ OptionKeys::fold_and_dock::slide_contact_frequency ]());

	random_mover.add_mover( move_anchor_mover, move_anchor_weight );
	random_mover.add_mover( rb_trial_mover, rb_weight );
	random_mover.add_mover( slide_mover, slide_weight );
}

bool AsymFoldandDockClaimer::read_tag( std::string tag, std::istream& is ) {
	if ( tag == "loop_file" || tag == "LOOP_FILE" ) {
		std::string file;
		is >> file;
		std::ifstream infile( file.c_str() );

		if ( !infile.good() ) {
			utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + file + "'" );
		}
		loops::PoseNumberedLoopFileReader loop_file_reader;
		loops::SerializedLoopList loops = loop_file_reader.read_pose_numbered_loops_file(infile, file, false );
		loops::Loops loop_defs = loops::Loops( loops ); // <==
		//  loop_defs = loop_defs.invert( input_pose_.total_residue() );
		tr << "Flexible residues: " << loop_defs << std::endl;
		moving_res_ = loop_defs;
	} else if ( tag == "CHAIN_BREAK_ASSYM_FND"  || tag == "chain_break_assym_fnd" ) {
		is >> chain_break_res_;
	} else  {
		throw utility::excn::EXCN_BadInput( " Unknown tag, only LOOP, CHAIN_BREAK_ASSYM_FND definition is allowed at this stage");
	}
	return true;
}

void AsymFoldandDockClaimer::initialize_dofs(
	core::pose::Pose& pose,
	claims::DofClaims const& init_dofs,
	claims::DofClaims& /*failed_to_init*/ ) {

	using namespace loops;
	using namespace kinematics;

	if ( moving_res_.size() == 0 ) throw utility::excn::EXCN_BadInput( " missing definition of moving residues, add a LOOP definition ");
	if ( moving_res_.size() > 1 ) throw utility::excn::EXCN_BadInput( " Only one movable region possible at this stage ");
	if ( chain_break_res_ == 0 ) throw utility::excn::EXCN_BadInput( " No chainbreak defined... ");

	/*
	Size moving_end( 1 );
	for ( Loops::const_iterator loop_it = moving_res_.begin(); loop_it != moving_res_.end(); ++loop_it ) {
	moving_start_ = loop_it->start();
	moving_end = loop_it->stop();
	}
	*/

	// make a PDBInfo for the pose...
	pose::PDBInfoOP pdb_info( new pose::PDBInfo( pose, true ) );

	utility::vector1< char > chain;
	for ( Size i=1; i<= pdb_info->nres(); ++i ) {
		if ( i <= chain_break_res_ ) {
			chain.push_back( 'A' );
		} else {
			chain.push_back( 'B' );
		}
	}

	pdb_info->set_chains( chain );
	pose.pdb_info( pdb_info );

	docking_local_refine_ = basic::options::option[ basic::options::OptionKeys::docking::docking_local_refine ]();
	bool slide ( !docking_local_refine_ );

	protocols::docking::DockingProtocol dock;
	protocols::docking::DockingInitialPerturbation dock_init(slide);
	std::string chainID("A_B");
	// If we don't have a docking jump we have to create it before calling setup_foldtree (strangely enough)
	if ( pose.fold_tree().num_jump() == 0 ) {
		FoldTree f(pose.fold_tree());
		f.clear();
		f.simple_tree( pose.total_residue() );
		f.new_jump( 1, pose.total_residue() , chain_break_res_ );
		pose.fold_tree( f );
	}
	protocols::docking::setup_foldtree( pose, chainID, dock.movable_jumps() );
	docking_jump_ = docking_jump( pose, chain_break_res_ );
	dock_init.apply( pose );
	// Setup the movemap
	kinematics::MoveMapOP movemap( new kinematics::MoveMap() );
	movemap->set_bb( true );
	movemap->set_jump( true );

	for ( claims::DofClaims::const_iterator it = init_dofs.begin(), eit = init_dofs.end();
			it != eit; ++it ) {
		if ( (*it)->owner().lock().get() == this ) {
			(*it)->toggle( *movemap, true );
		}
	}

}

void AsymFoldandDockClaimer::generate_claims( claims::DofClaims& new_claims ) {
	// Set all cuts to real cuts. We don't want to close any of them...
	loops::Loops loops;
	loops::Loop single_loop;

	std::cout << input_pose_.fold_tree() << std::endl;
	utility::vector1< int > cuts( input_pose_.conformation().fold_tree().cutpoints() );
	for ( Size i = 1; i <= cuts.size(); ++i ) {

		new_claims.push_back( claims::DofClaimOP( new claims::CutClaim( get_self_weak_ptr(), std::make_pair( Parent::label(), cuts[i]),
			claims::DofClaim::INIT // for now... eventually CAN_INIT ?
			) ) );
	}
}

core::Size AsymFoldandDockClaimer::docking_jump( core::pose::Pose& pose, core::Size chain_break_res )
{
	using namespace core::kinematics;

	core::Size nres ( pose.total_residue() );
	// Does the chain start at a reasonable place in the sequence?
	runtime_assert( chain_break_res > 1 && chain_break_res <= nres );

	// find the anchor
	FoldTree f( pose.fold_tree() );
	int jump_number(1);
	for ( Size i=1; i<= f.num_jump(); ++i ) {
		Size res( f.downstream_jump_residue( i ) );
		if ( res > chain_break_res ) {
			jump_number = i;
		}
	}
	return jump_number;
}

} //topology_broker
} //protocols
