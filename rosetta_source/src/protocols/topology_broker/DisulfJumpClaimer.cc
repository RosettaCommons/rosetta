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
/// @detailed responsibilities:
/// @author Robert Vernon

// Unit Headers
#include <protocols/topology_broker/DisulfJumpClaimer.hh>

// Package Headers
#include <protocols/topology_broker/DofClaim.hh>
#include <protocols/jumping/DisulfPairingsList.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
// AUTO-REMOVED #include <core/fragment/FragID.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/FrameList.hh>
// AUTO-REMOVED #include <core/fragment/FragmentIO.hh>
#include <protocols/basic_moves/FragmentMover.hh>
// AUTO-REMOVED #include <protocols/jumping/JumpSetup.hh>
#include <core/fragment/JumpingFrame.hh>
// AUTO-REMOVED #include <core/fragment/BBTorsionSRFD.hh>
// AUTO-REMOVED #include <core/fragment/JumpSRFD.hh>

#include <core/fragment/FragData.fwd.hh>



// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


//#include <basic/options/option.hh>

//// C++ headers

// option key includes

static basic::Tracer tr("protocols.topo_broker",basic::t_info);
//static numeric::random::RandomGenerator RG(18828234);

namespace protocols {
namespace topology_broker {

using namespace core;

DisulfJumpClaimer::DisulfJumpClaimer() :
	//jump_def_( NULL ),
	//	init_mover_( NULL ),
	bKeepJumpsFromInputPose_( true )
{
	set_bInitDofs( true ); //we want to initialize jumps
}

// DisulfJumpClaimer::DisulfJumpClaimer( std::string const& tag, weights::AbinitioMoverWeightOP weight ) :
// 	FragmentClaimer( NULL, tag, weight ),
// 	//	jump_def_ ( jump_def ),
// 	//	init_mover_( NULL ),
// 	bKeepJumpsFromInputPose_( true )
// {
// 	set_bInitDofs( true ); //we want to initialize jumps
// }

void DisulfJumpClaimer::new_decoy() {
	//if ( !jump_def_ ) return;
  //runtime_assert( jump_def_ );

	generate_jump_frags( *jumping::StandardDisulfPairingLibrary::get_instance(), all_frames_ );

 	core::fragment::FragSetOP jump_frags = new core::fragment::OrderedFragSet;
 	jump_frags->add( all_frames_ );

	basic_moves::ClassicFragmentMoverOP mover;
  mover = new basic_moves::ClassicFragmentMover( jump_frags, movemap_ );
	mover->type( mover_tag() );
  mover->set_check_ss( false ); // this doesn't make sense with jump fragments
	mover->enable_end_bias_check( false ); //no sense for discontinuous fragments
	set_mover( mover );

	//Size attempts( 10 );
	//do {
	//	current_jumps_ = jump_def_->create_jump_sample();
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

void DisulfJumpClaimer::initialize_dofs( core::pose::Pose& pose, DofClaims const& init_dofs, DofClaims& failed_to_init ) {

	//init_mover_ = new basic_moves::ClassicFragmentMover( jump_frags, movemap_ );
	//init_mover_->type( mover_tag() );
 	//init_mover_->set_check_ss( false ); // this doesn't make sense with jump fragments
 	//init_mover_->enable_end_bias_check( false ); //no sense for discontinuous fragments

// 	kinematics::MoveMapOP movemap = new kinematics::MoveMap();

// 	for ( DofClaims::const_iterator it = init_dofs.begin(), eit = init_dofs.end();
// 				it != eit; ++it ) {
//     if ( (*it)->owner()==this ) {
//       (*it)->toggle( *movemap, true );
// 		}
// 	}




	//need to copy coords and jumps --- if chunks were idealized no problem .... but non-idealized stuff ?
// 	if ( init_mover_ ) {
// 		basic_moves::FragmentMoverOP frag_mover = get_frag_mover_ptr();
// 		set_mover( init_mover_ );
// 		FragmentClaimer::initialize_dofs( pose, init_dofs, failed_to_init );
// 		set_mover( frag_mover );
// 		init_mover_ = NULL;
// 	} else {
 		FragmentClaimer::initialize_dofs( pose, init_dofs, failed_to_init );
// 	}
}

void DisulfJumpClaimer::generate_jump_frags(
	protocols::jumping::DisulfPairingLibrary const& lib,
	//core::kinematics::MoveMap const& mm,
	core::fragment::FrameList& all_frames
) const
{
	typedef utility::vector1< Size > JumpList;
	all_frames.reserve( all_frames.size() + all_jump_pairings_.size() );

	core::fragment::FragDataList frag_data;
	lib.create_jump_fragments( false, frag_data );

	for ( Size i = 1; i <= all_jump_pairings_.size(); ++i) {
		int const startpos( all_jump_pairings_[i].pos1 );
		int const endpos( all_jump_pairings_[i].pos2 );

		core::fragment::JumpingFrameOP frame =
			new core::fragment::JumpingFrame( startpos, endpos, 2 );

		frame->set_pos( 1, startpos );
		frame->set_pos( 2, endpos );

		frame->add_fragment( frag_data );

		runtime_assert( frame->nr_frags() );
		all_frames.push_back( frame );
	}

}

// void DisulfJumpClaimer::generate_jump_frames(
// 	 core::fragment::FrameList& all_frames,
// 	 core::kinematics::MoveMap const& mm
// ) const
// {
// 	all_frames.reserve( all_frames.size() + all_jump_pairings_.size() );

// 	for ( Size i = 1; i <= all_jump_pairings_.size(); ++i) {
// 		int const startpos( all_jump_pairings_[i].pos1 );
// 		int const endpos( all_jump_pairings_[i].pos2 );

// 		core::fragment::FragDataOP frag_data = new core::fragment::FragData;

// 		frag_data->add_residue( new core::fragment::UpJumpSRFD );
// 		frag_data->add_residue( new core::fragment::DownJumpSRFD );

// 		core::fragment::JumpingFrameOP frame =
// 			new core::fragment::JumpingFrame( startpos, endpos, frag_data->size () );

// 		Size pos = 1;
// 		frame->set_pos( pos++, startpos );
// 		frame->set_pos( pos++, endpos );

// 		frame->add_fragment( frag_data );

// 		runtime_assert( frame->nr_frags() );
// 		all_frames.push_back( frame );
// 	}
// }


void DisulfJumpClaimer::generate_claims( DofClaims& new_claims ) {

	// get flexible jumps ( beta-sheet stuff etc. )
	/// in future get rid of JumpSample class all-together.

	movemap_->set_jump( true ); //we switch them off on a as-need basis
	movemap_->set_bb( true );

// 	core::fragment::FragSetOP jump_frags = new core::fragment::OrderedFragSet;
// 	core::fragment::FrameList jump_frames;

// 	//generate_jump_frames( jump_frames, *movemap_ );
// 	jump_frags->add( jump_frames );


	for ( Size i = 1; i <= all_jump_pairings_.size(); ++i) {

		Size const up( all_jump_pairings_[i].pos1 );
		Size const down( all_jump_pairings_[i].pos2 );

			//new_claims.push_back( new JumpClaim( this, up, down, up_atom, down_atom, DofClaim::INIT ) );
		new_claims.push_back( new JumpClaim( this, up, down, DofClaim::INIT ) );

	}

}

bool DisulfJumpClaimer::read_tag( std::string tag, std::istream& is ) {
	if ( tag == "NO_USE_INPUT_POSE" ) {
		bKeepJumpsFromInputPose_ = false;
	} else if ( tag == "DISULF" ) {
		Size pos1, pos2;
		char ss1, ss2;

		is >> pos1 >> ss1 >> pos2 >> ss2;

		protocols::jumping::DisulfPairing dis_pair;

		dis_pair.pos1 = pos1;
		dis_pair.pos2 = pos2;
		dis_pair.seq_sep = abs((int)(pos1-pos2));
		dis_pair.ss_type = 1;//PLACEHOLDER!

		all_jump_pairings_.push_back( dis_pair );

	} else if ( tag == "mover_weight" ) {
		read_mover_weight( is );
	} else return Parent::read_tag( tag, is );
	return true;
}

} //topology_broker
} //protocols
