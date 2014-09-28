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
/// @detailed responsibilities:
/// @author Oliver Lange

// Unit Headers
#include <protocols/topology_broker/DensityScoringClaimer.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/BBClaim.hh>
#include <protocols/topology_broker/claims/LegacyRootClaim.hh>
#include <protocols/topology_broker/claims/JumpClaim.hh>
#include <protocols/topology_broker/Exceptions.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/SequenceNumberResolver.hh>


// Project Headers
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>

// AUTO-REMOVED #include <protocols/jumping/ResiduePairJumpSetup.hh>
// AUTO-REMOVED #include <protocols/jumping/ResiduePairJump.hh>
// AUTO-REMOVED #include <protocols/jumping/JumpSetup.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
// AUTO-REMOVED #include <core/sequence/util.hh>


#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
//#include <basic/options/option.hh>

//// C++ headers
// AUTO-REMOVED #include <fstream>

#include <protocols/jumping/JumpSetup.fwd.hh>
#include <utility/vector1.hh>


// option key includes


static thread_local basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace core;
DensityScoringClaimer::DensityScoringClaimer()
{}

/* void DensityScoringClaimer::initialize_residues( core::pose::Pose& pose, claims::SequenceClaimOP my_claim, claims::DofClaims& failed_to_init ) {
	SequenceClaimer::initialize_residues( pose, my_claim, failed_to_init );
	resolved_anchor_residue_ = broker().sequence_number_resolver().find_global_pose_number( anchor_chain_, anchor_residue_ );
	//resolved_anchor_residue_ = broker().resolve_residue( anchor_chain_, anchor_residue_ );
	vrt_id_ = broker().sequence_number_resolver().find_global_pose_number(my_claim->label());
	//vrt_id_ = my_claim->offset(); // now that we know it, set vrt id here
	tr.Debug << "Setting vrt_id_ to " << vrt_id_ << std::endl;
	tr.Debug << "Setting resolved_anchor_residue_ to " << resolved_anchor_residue_ << std::endl;
} */


void DensityScoringClaimer::generate_sequence_claims( claims::DofClaims& new_claims ){
	//SequenceClaim for ElectronDensityCenter
	new_claims.push_back ( claims::DofClaimOP( new claims::SequenceClaim( get_self_weak_ptr(), "X", label() ) ));
}


void DensityScoringClaimer::generate_claims( claims::DofClaims& new_claims ) {

	std::pair < std::string, core::Size> local_vrt_pos ( label() ,1);
	std::pair < std::string, core::Size> local_anchor_pos ( anchor_chain_, anchor_residue_ );
	new_claims.push_back( claims::DofClaimOP( new claims::LegacyRootClaim( get_self_weak_ptr(), local_vrt_pos, claims::DofClaim::EXCLUSIVE) ) );
	new_claims.push_back( claims::DofClaimOP( new claims::JumpClaim( get_self_weak_ptr(), local_vrt_pos, local_anchor_pos, "ORIG", "CA", claims::DofClaim::EXCLUSIVE ) ) );

	//Claim the BB torsion for the virtual residue, which shouldn't be modified anyway, to make the Broker happy
	new_claims.push_back( claims::DofClaimOP( new claims::BBClaim( get_self_weak_ptr(), std::make_pair( label(), 1 ) ) ) );

}

TopologyClaimerOP DensityScoringClaimer::clone() const {
	return TopologyClaimerOP( new DensityScoringClaimer( *this ) );
}


void DensityScoringClaimer::add_constraints( core::pose::Pose& ) const {
	// ??
}

void  DensityScoringClaimer::set_defaults() {
	set_label( "ElectronDensityCenter" );	// Needed for SequenceClaim generated with generate_sequence_claims for Residue X.
	anchor_chain_ = ""; //usually anchored to DEFAULT chain
	anchor_residue_ = 0;
}

bool DensityScoringClaimer::read_tag( std::string tag, std::istream& is ) {
	using namespace jumping;
	if ( tag == "anchor_residue" ) {
		is >> anchor_residue_;
		tr.Debug << "Read anchor res = " << anchor_residue_ << std::endl;
	} else if ( tag == "anchor_chain" ){
		is >> anchor_chain_;
		tr.Debug << "Reading anchor chain: " << anchor_chain_ << std::endl;
	}
	/*else if ( SequenceClaimer::read_tag( tag, is ) ) {
		//noop
	} */ else return false;
	return true;
}

bool DensityScoringClaimer::accept_declined_claim( claims::DofClaim const& was_declined ) {
	tr.Error << "JumpClaim of " << /*<< type() <<*/ " was declined: " << was_declined << std::endl;
	return false; // full tolerance here ---
}


void DensityScoringClaimer::init_after_reading() {
	if ( anchor_residue_ == 0 || anchor_chain_ == "" ) {
		throw EXCN_Input( "need to specify both anchor residue (tag = anchor_residue) and anchor chain (tag = anchor_chain) for DensityScoring" );
	}

	/*set_label( ObjexxFCL::string_of( anchor_residue_ ) + anchor_chain_ );
	set_sequence( "X" );
	set_label( "ElectronDensityCenter" );
	*/
}



} //topology_broker
} //protocols
