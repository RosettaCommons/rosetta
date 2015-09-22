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
#include <protocols/topology_broker/BasicJumpClaimer.hh>

// Project Headers
#include <protocols/topology_broker/claims/SequenceClaim.hh>
#include <protocols/topology_broker/SequenceNumberResolver.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/FoldTree.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

//C++ Headers
#include <string>
#include <sstream>
#include <utility>

static THREAD_LOCAL basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

void BasicJumpClaimer::generate_claims( claims::DofClaims& claims ){
	//Verify that jump is possible
	claims::SequenceClaim start_chain_claim = broker().resolve_sequence_label( start_label_ );
	claims::SequenceClaim end_chain_claim   = broker().resolve_sequence_label( end_label_   );

	if ( ! ( ( start_position_ > 0 ) && ( start_position_ <= start_chain_claim.length() ) ) ) {
		std::ostringstream msg;
		msg << "Jump start position " << start_position_ << " in label '" << start_label_
			<< "' requested in BasicJumpClaimer is not a valid position within that label. Bounds are ["
			<< broker().sequence_number_resolver().find_global_pose_number( start_label_ ) << ", "
			<< broker().sequence_number_resolver().find_global_pose_number( start_label_ , start_chain_claim.length() )
			<< "]." << std::endl;
		throw utility::excn::EXCN_BadInput( msg.str() );
	} else if ( ! ( ( end_position_ > 0 ) && ( end_position_ <= end_chain_claim.length() ) ) ) {
		std::ostringstream msg;
		msg << "Jump end position " << end_position_ << " in label '" << end_label_
			<< "' requested in BasicJumpClaimer is not a valid position within that label. Bounds are ["
			<< broker().sequence_number_resolver().find_global_pose_number( end_label_ ) << ", "
			<< broker().sequence_number_resolver().find_global_pose_number( end_label_ , end_chain_claim.length() )
			<< "]." << std::endl;
		throw utility::excn::EXCN_BadInput( msg.str() );
	}

	//Add new jump claim
	if ( ( start_atom_ == "" ) && ( end_atom_ == "" ) ) {
		claims.push_back( claims::DofClaimOP( new claims::JumpClaim( get_self_weak_ptr(), std::make_pair( start_label_, start_position_ ),
			std::make_pair( end_label_, end_position_ ) ) ) );
	} else if ( ( start_atom_ != "" ) && ( end_atom_ != "" ) ) {
		claims.push_back( claims::DofClaimOP( new claims::JumpClaim( get_self_weak_ptr(), std::make_pair( start_label_, start_position_ ),
			std::make_pair( end_label_, end_position_ ), start_atom_,
			end_atom_ ) ) );
	} else {
		std::ostringstream msg;
		msg << "BasicJumpClaimer with JumpClaim between '" << start_label_ << "' (position "
			<< start_position_ << ") and " << end_label_ << " (position " << end_position_
			<< " did not set both (or neither) jump atoms." << std::endl;
		throw utility::excn::EXCN_BadInput( msg.str() );
	}
}

void BasicJumpClaimer::initialize_dofs( core::pose::Pose& pose, claims::DofClaims const& claims, claims::DofClaims&)
{
	claims::JumpClaimOP our_claim;
	for ( claims::DofClaims::const_iterator claim = claims.begin(); claim != claims.end(); ++claim ) {
		if ( claim->get()->owner().lock().get() == this ) {
			our_claim = utility::pointer::dynamic_pointer_cast< claims::JumpClaim >( *claim );
			runtime_assert( our_claim != 0 );
		}
	}

	core::Size jump_number = pose.fold_tree().
		jump_nr( broker().sequence_number_resolver().find_global_pose_number( our_claim->local_pos1() ),
		broker().sequence_number_resolver().find_global_pose_number( our_claim->local_pos2()   ) );
	runtime_assert( jump_number );

	core::kinematics::Jump init_jump;
	//Set jump (initialized to zero by default constructor) to length 50 in a random direction
	init_jump.random_trans( 50 );

	pose.set_jump( (int) jump_number, init_jump );
}

bool BasicJumpClaimer::read_tag( std::string tag, std::istream& is ){
	// Expects input JUMP_[START,END] label internal_position_number atom_letter
	// The atom letter is optional

	std::string line;
	if ( ( tag == "JUMP_START" ) || ( tag == "jump_start" ) || ( tag == "Jump_Start" ) ) {
		getline( is, line );
		tr.Debug << "BasicJumpClaimer read line " << line << std::endl;
		std::stringstream linestringstream( line );

		linestringstream >> start_label_;
		linestringstream >> start_position_;
		if ( ! ( linestringstream >> start_atom_ ) ) {
			start_atom_ = "";
		}
	} else if  ( ( tag == "JUMP_END" ) || ( tag == "jump_end" ) || ( tag == "Jump_End") ) {
		getline( is, line );
		tr.Debug << "BasicJumpClaimer read line " << line << std::endl;
		std::stringstream linestringstream( line );

		linestringstream >> end_label_;
		linestringstream >> end_position_;
		if ( ! ( linestringstream >> end_atom_ ) ) {
			end_atom_ = "";
		}
	} else {
		tr.Debug << "BasicJumpClaimer calling parent function on tag " << tag << std::endl;
		return Parent::read_tag( tag, is );
	}

	return true;
}

}
}
