// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DisulfJumpClaimer
/// @brief  Claimer for disulfide jump sampling
/// @details responsibilities:
/// @author Robert Vernon

// Unit Headers
#include <protocols/topology_broker/DisulfJumpClaimer.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/JumpClaim.hh>
#include <protocols/jumping/DisulfPairingsList.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/FrameList.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <core/fragment/JumpingFrame.hh>
#include <core/fragment/FragData.fwd.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

static thread_local basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace core;

DisulfJumpClaimer::DisulfJumpClaimer() :
	//jump_def_( NULL ),
	// init_mover_( NULL ),
	bKeepJumpsFromInputPose_( true )
{
	set_bInitDofs( true ); //we want to initialize jumps
}

DisulfJumpClaimer::~DisulfJumpClaimer() {}

TopologyClaimerOP DisulfJumpClaimer::clone() const {
	return TopologyClaimerOP( new DisulfJumpClaimer( *this ) );
}

void DisulfJumpClaimer::new_decoy() {
	generate_jump_frags( *jumping::StandardDisulfPairingLibrary::get_instance(), all_frames_ );

	core::fragment::FragSetOP jump_frags( new core::fragment::OrderedFragSet );
	jump_frags->add( all_frames_ );

	simple_moves::ClassicFragmentMoverOP mover;
	mover = simple_moves::ClassicFragmentMoverOP( new simple_moves::ClassicFragmentMover( jump_frags, movemap_ ) );
	mover->type( mover_tag() );
	mover->set_check_ss( false ); // this doesn't make sense with jump fragments
	mover->enable_end_bias_check( false ); //no sense for discontinuous fragments
	set_mover( mover );

	//Size attempts( 10 );
	//do {
	// current_jumps_ = jump_def_->create_jump_sample();
	//} while ( !current_jumps_.is_valid() && attempts-- );

	//if ( !current_jumps_.is_valid() ) {
	//    utility_exit_with_message( "not able to build valid fold-tree in DisulfJumpClaimer" );
	//}
	//tr.Debug << "current_jumps " << current_jumps_ << std::endl;
}


void DisulfJumpClaimer::new_decoy( core::pose::Pose const& pose ) {
	core::pose::Pose temp(pose);

	//Figure out what to do when the pose already has a fold tree later
	new_decoy();
}

void DisulfJumpClaimer::initialize_dofs( core::pose::Pose& pose, claims::DofClaims const& init_dofs, claims::DofClaims& failed_to_init ) {

	//init_mover_ = new simple_moves::ClassicFragmentMover( jump_frags, movemap_ );
	//init_mover_->type( mover_tag() );
	//init_mover_->set_check_ss( false ); // this doesn't make sense with jump fragments
	//init_mover_->enable_end_bias_check( false ); //no sense for discontinuous fragments

	//  kinematics::MoveMapOP movemap = new kinematics::MoveMap();

	//  for ( DofClaims::const_iterator it = init_dofs.begin(), eit = init_dofs.end();
	//     it != eit; ++it ) {
	//     if ( (*it)->owner()==this ) {
	//       (*it)->toggle( *movemap, true );
	//   }
	//  }


	//need to copy coords and jumps --- if chunks were idealized no problem .... but non-idealized stuff ?
	//  if ( init_mover_ ) {
	//   simple_moves::FragmentMoverOP frag_mover = get_frag_mover_ptr();
	//   set_mover( init_mover_ );
	//   FragmentClaimer::initialize_dofs( pose, init_dofs, failed_to_init );
	//   set_mover( frag_mover );
	//   init_mover_ = NULL;
	//  } else {
	FragmentClaimer::initialize_dofs( pose, init_dofs, failed_to_init );
	//  }
}

void DisulfJumpClaimer::generate_jump_frags(
	protocols::jumping::DisulfPairingLibrary const& lib,
	//core::kinematics::MoveMap const& mm,
	core::fragment::FrameList& all_frames
) const
{
	all_frames.reserve( all_frames.size() + local_disulf_data_.size() );

	//TODO False assert. Fragments have to be dealt with properly here and idk how to do it atm.
	runtime_assert( false );
	core::fragment::FragDataOPs frag_data;
	lib.create_jump_fragments( false, frag_data );

	for ( Size i = 1; i <= all_jump_pairings_.size(); ++i ) {
		int const startpos( all_jump_pairings_[i].pos1 );
		int const endpos( all_jump_pairings_[i].pos2 );

		core::fragment::JumpingFrameOP frame( new core::fragment::JumpingFrame( startpos, endpos, 2 ) );

		frame->set_pos( 1, startpos );
		frame->set_pos( 2, endpos );

		frame->add_fragment( frag_data );

		runtime_assert( frame->nr_frags() );
		all_frames.push_back( frame );
	}

}

// void DisulfJumpClaimer::generate_jump_frames(
//   core::fragment::FrameList& all_frames,
//   core::kinematics::MoveMap const& mm
// ) const
// {
//  all_frames.reserve( all_frames.size() + all_jump_pairings_.size() );

//  for ( Size i = 1; i <= all_jump_pairings_.size(); ++i) {
//   int const startpos( all_jump_pairings_[i].pos1 );
//   int const endpos( all_jump_pairings_[i].pos2 );

//   core::fragment::FragDataOP frag_data = new core::fragment::FragData;

//   frag_data->add_residue( new core::fragment::UpJumpSRFD );
//   frag_data->add_residue( new core::fragment::DownJumpSRFD );

//   core::fragment::JumpingFrameOP frame =
//    new core::fragment::JumpingFrame( startpos, endpos, frag_data->size () );

//   Size pos = 1;
//   frame->set_pos( pos++, startpos );
//   frame->set_pos( pos++, endpos );

//   frame->add_fragment( frag_data );

//   runtime_assert( frame->nr_frags() );
//   all_frames.push_back( frame );
//  }
// }


void DisulfJumpClaimer::generate_claims( claims::DofClaims& new_claims ) {
	using claims::LocalPosition;
	using std::pair;

	//Initialize all_jump_pairings_ list with the data gathered during read_tag
	for ( utility::vector1< claims::JumpClaimOP >::iterator bond_it = local_disulf_data_.begin();
			bond_it != local_disulf_data_.end(); ++bond_it ) {
		claims::JumpClaimOP claim( *bond_it );

		core::Size pos1 = claim->global_pos1();
		core::Size pos2 = claim->global_pos2();

		protocols::jumping::DisulfPairing dis_pair;

		dis_pair.pos1 = pos1;
		dis_pair.pos2 = pos2;
		dis_pair.seq_sep = std::abs((int)(pos1-pos2));
		//TODO: use ss1, ss2 to do... something? with the secondary structure
		dis_pair.ss_type = 1;//PLACEHOLDER!

		all_jump_pairings_.push_back( dis_pair );

		new_claims.push_back( claims::DofClaimOP( new claims::JumpClaim( get_self_weak_ptr(),
			claim->local_pos1(),
			claim->local_pos2(),
			claims::DofClaim::INIT ) ) );
	}

	// get flexible jumps ( beta-sheet stuff etc. )
	/// in future get rid of JumpSample class all-together.

	movemap_->set_jump( true ); //we switch them off on a as-need basis
	movemap_->set_bb( true );

	//  core::fragment::FragSetOP jump_frags = new core::fragment::OrderedFragSet;
	//  core::fragment::FrameList jump_frames;

	//  //generate_jump_frames( jump_frames, *movemap_ );
	//  jump_frags->add( jump_frames );
}

bool DisulfJumpClaimer::read_tag( std::string tag, std::istream& is ) {
	if ( tag == "NO_USE_INPUT_POSE" ) {
		bKeepJumpsFromInputPose_ = false;
	} else if ( tag == "DISULF" ) {
		//Requires syntax label1 pos1 ss1 label2 pos2 ss2
		std::string label1, label2;
		Size pos1, pos2;
		std::string ss1, ss2;

		is >> pos1 >> label1 >> ss1 >> pos2 >> label2 >> ss2;

		if ( !( ss1 == "S" || ss1 == "H" || ss1 == "E" )  ) {
			throw utility::excn::EXCN_BadInput(
				"When reading DisulfJumpClaimer, secondary structure character '"
				+ss1+"' was invalid. Valid characters are 'S', 'H', and 'E'." );
		} else if ( !( ss2 == "S" || ss2 == "H" || ss2 == "E" )  ) {
			throw utility::excn::EXCN_BadInput(
				"When reading DisulfJumpClaimer, secondary structure character '"
				+ss2+"' was invalid. Valid characters are 'S', 'H', and 'E'." );
		}

		claims::LocalPosition local_pos1 = std::make_pair( label1, pos1 );
		claims::LocalPosition local_pos2 = std::make_pair( label2, pos2 );

		//Use jump claim's atom to keep track of the secondary structure (a bit hacky, I know)
		claims::JumpClaimOP disulf_bond( new claims::JumpClaim( get_self_weak_ptr(),
			std::make_pair( label1, pos1 ),
			std::make_pair( label2, pos2 ),
			ss1,
			ss2 ) );

		local_disulf_data_.push_back( disulf_bond );

	} else if ( tag == "mover_weight" ) {
		read_mover_weight( is );
	} else return Parent::read_tag( tag, is );
	return true;
}

} //topology_broker
} //protocols
