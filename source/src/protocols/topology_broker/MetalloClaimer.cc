// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange

// Unit Headers
#include <protocols/topology_broker/MetalloClaimer.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/SequenceClaim.hh>
#include <protocols/topology_broker/Exceptions.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/SequenceNumberResolver.hh>


// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/Pose.hh>

#include <protocols/jumping/ResiduePairJumpSetup.hh>
#include <protocols/jumping/ResiduePairJump.hh>
#include <protocols/jumping/JumpSetup.hh>

// ObjexxFCL Headers

// Utility headers


//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
//#include <basic/options/option.hh>

//// C++ headers

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace core;
MetalloClaimer::MetalloClaimer() : residue_pair_jump_( /* NULL */ )
{}

//bool MetalloClaimer::allow_claim( DofClaim const& foreign_claim ) {
// 	if ( foreign_claim.owner() == this ) return true; // always allow your own claims!
// 	if ( foreign_claim.type() == DofClaim::SEQUENCE ) {
// 		if ( foreign_claim.pos( 1 ) == 1  && foreing_claim.right() == DofClaim::EXCLUSIVE ) {
// 			return false; // don't accept any other fixed positions right now --- maybe never ?!
// 		}
// 	}
//	return true;
//} // MetalloClaimer::allow_claim()

//void MetalloClaimer::initialize_residues( core::pose::Pose& pose, claims::SequenceClaimOP my_claim, claims::DofClaims& failed_to_init ) {
void MetalloClaimer::generate_claims( claims::DofClaims& new_claims ) {

	// TODO: refactor FragmentJumpClaimer to deal with this correctly.
	// What we want to do is create a jump from the zinc to the anchor chain, and build the appropriate fragments including
	// the jump itself? This is done with old "jump interval" code that is nasty nasty.
	// Also, much of this code is broken/stupid, because we're in the process of refactoring.
	//runtime_assert( false );

	//Find claim of anchor chain
	claims::SequenceClaim my_claim = broker().resolve_sequence_label( label() );

	jump_setup_->clear();
	resolved_anchor_residue_ = broker().sequence_number_resolver().find_global_pose_number( anchor_chain_, anchor_residue_ );
	core::Size my_claim_pos = broker().sequence_number_resolver().find_global_pose_number(my_claim.label());

	tr.Trace << "MetalloClaimer: setup jump between " << resolved_anchor_residue_ << " " << my_claim_pos<< std::endl;

	core::Size cut_interval_low  = my_claim_pos - 1;
	core::Size cut_interval_high = my_claim_pos;

	std::pair<core::Size, std::string> terminal_label = broker().sequence_number_resolver().terminal_pair();
	core::Size seq_length = terminal_label.first + broker().resolve_sequence_label( terminal_label.second ).length();

	if (cut_interval_low < 1) {
		cut_interval_low = 1; // no cut in front of first element of sequence possible
	}
	if (cut_interval_high >= seq_length ) {
		cut_interval_high = seq_length -1;	// no cut possible behind the last residue
	}


	tr.Trace << "Cut Interval: " << cut_interval_low << " " << cut_interval_high << std::endl;

	//Size cut_position = my_claim->offset() +

	core::Size jump_start;
	core::Size jump_end;

	if ( resolved_anchor_residue_ < my_claim_pos ) {
		jump_start = resolved_anchor_residue_;
		jump_end = my_claim_pos;
	} else {
		jump_start = my_claim_pos;
		jump_end = resolved_anchor_residue_;
	}

	jump_setup_->add_jump(
			jumping::Interval( jump_start, jump_end ), //jump
			jumping::Interval( cut_interval_low, cut_interval_high ) //cutpoint interval
	);
	new_decoy();

	FragmentJumpClaimer::generate_claims( new_claims,  anchor_chain_, label() );
	SequenceClaimer::generate_claims( new_claims );
}

void MetalloClaimer::add_constraints( core::pose::Pose& ) const {
}

void  MetalloClaimer::set_defaults() {
	SequenceClaimer::set_defaults();
	FragmentJumpClaimer::set_defaults();
	anchor_chain_ = ""; //usually anchored to DEFAULT chain
	anchor_residue_ = 0;
	residue_pair_jump_ = jumping::ResiduePairJumpOP( new jumping::ResiduePairJump );
}

bool MetalloClaimer::read_tag( std::string tag, std::istream& is ) {
	using namespace jumping;
	if ( tag == "ligand" ) {
		is >> ligand_;
	} else if ( tag == "anchor" ) {
		is >> anchor_residue_;
	} else if ( tag == "aa" ) {
		for ( int i = 1; i <= 2; ++i ) {
			std::string name;
			core::chemical::ResidueTypeSetCOP residue_type_set(
					core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
			is >> name;
			if ( residue_type_set->name_map(name).is_protein() )
				name = name + ":CtermProteinFull:NtermProteinFull";
			core::chemical::ResidueType const & res_type( residue_type_set->name_map(name) );
			residue_pair_jump_->add_residue_single( res_type );
		}
	} else if ( tag == "cst_atoms" ) {
		for ( int i = 1; i <= 2; ++i ) {
			std::string name;
			for ( int j = 1; j <= 3; ++j ) {
				is >> name;
				residue_pair_jump_->set_cstAtoms(i,j,name);
			}
		}
	} else if ( tag == "jump_atoms" ) {
		for ( int i = 1; i <= 2; ++i ) {
			std::string name;
			for ( int j = 1; j <= 3; ++j ) {
				is >> name;
				residue_pair_jump_->set_jumpAtoms(i,j,name);
			}
		}
	} else if ( tag == "disAB"
			|| tag == "angleA" || tag == "angleB"
			|| tag == "dihedralA" || tag == "dihedralB" || tag == "dihedralAB" ) {
		//these are all multi valued read until line ends
		std::string line;
		getline( is, line );
		std::istringstream in( line );
		Real value;
		if ( tag == "disAB" ) {
			while (in >> value) {
				residue_pair_jump_->set_cstInfo( disAB, value );
			}
		} else if ( tag == "angleA" ) {
			while (in >> value) {
				residue_pair_jump_->set_cstInfo( angleA, value );
			}
		} else if ( tag == "angleB" ) {
			while (in >> value) {
				residue_pair_jump_->set_cstInfo( angleB, value );
			}
		} else if ( tag == "dihedralA" ) {
			while (in >> value) {
				residue_pair_jump_->set_cstInfo( dihedralA, value );
			}
		} else if ( tag == "dihedralB" ) {
			while (in >> value) {
				residue_pair_jump_->set_cstInfo( dihedralB, value );
			}
		} else if ( tag == "dihedralAB" ) {
			while (in >> value) {
				residue_pair_jump_->set_cstInfo( dihedralAB, value );
			}
		} else runtime_assert( 0 ); //if you are here you missed a tag in the statement above.
	} else if ( SequenceClaimer::read_tag( tag, is ) ) {
		//noop
	} else if ( FragmentJumpClaimer::read_tag( tag, is ) ) {
		//noop
	} else return false;
	return true;
}

void MetalloClaimer::init_after_reading() {
	if ( !anchor_residue_ ) {
		throw EXCN_Input( "need to specify anchor residue for MetalloLigand "+ligand_ );
	}

	if (anchor_chain_ == ""){
		tr.Info << "anchor chain has not been specified. DEFAULT chain will be used." << std::endl;
		anchor_chain_ = "DEFAULT";
	}

	set_label( ligand_ + ObjexxFCL::string_of( anchor_residue_ ) + anchor_chain_ );

	//Parent SequenceClaimer should only produce one claim, that label is the anchor_chain_.
	//runtime_assert( get_sequence_labels().size() == 1);
	//if ( !anchor_chain_.size() ) anchor_chain_ = get_sequence_labels().at(1);
	set_sequence( "Z["+ligand_+"]" );
	//set_label( "metal_"+ligand_ );

	//start_label( anchor_chain_ );
	//start_position( anchor_residue_ );
	//end_label(label() );
	//end_position ( 1 );

	residue_pair_jump_->init_mini_pose();
	jump_setup_ = jumping::ResiduePairJumpSetupOP( new jumping::ResiduePairJumpSetup );
	jump_setup_->add_residue_pair_jump( residue_pair_jump_ );
	set_jump_def( jump_setup_ ); //tell the underlying FragmentJumpClaimer
}


} //topology_broker
} //protocols
