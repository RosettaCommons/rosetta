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
#include <protocols/topology_broker/DofClaim.hh>
#include <protocols/topology_broker/Exceptions.hh>
#include <protocols/topology_broker/TopologyBroker.hh>

// Project Headers
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED
#include <core/pose/Pose.hh>
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

//Auto Headers
#include <protocols/jumping/JumpSetup.fwd.hh>


// option key includes


static basic::Tracer tr("protocols.topo_broker",basic::t_info);

namespace protocols {
namespace topology_broker {

using namespace core;
DensityScoringClaimer::DensityScoringClaimer()
{}

void DensityScoringClaimer::initialize_residues( core::pose::Pose& pose, SequenceClaimOP my_claim, DofClaims& failed_to_init ) {
	SequenceClaimer::initialize_residues( pose, my_claim, failed_to_init );
	resolved_anchor_residue_ = broker().resolve_residue( anchor_chain_, anchor_residue_ );

	vrt_id_ = my_claim->offset(); // now that we know it, set vrt id here
	tr.Debug << "Setting vrt_id_ to " << vrt_id_ << std::endl;
	tr.Debug << "Setting resolved_anchor_residue_ to " << resolved_anchor_residue_ << std::endl;
}

void DensityScoringClaimer::generate_claims( DofClaims& new_claims ) {
	new_claims.push_back( new RootClaim( this, vrt_id_, DofClaim::EXCLUSIVE ) );
	new_claims.push_back( new JumpClaim( this, vrt_id_, resolved_anchor_residue_, "ORIG", "CA", DofClaim::EXCLUSIVE ) );
	//new_claims.push_back( new JumpClaim( this, resolved_anchor_residue_, vrt_id_,  "ORIG", "CA", DofClaim::EXCLUSIVE ) );

	// sequence claimer does not claim the cut automatically??
	//   ??? it should ....
	//new_claims.push_back( new CutClaim( this, vrt_id_-1, DofClaim::EXCLUSIVE ) );

	SequenceClaimer::generate_claims( new_claims );
}

TopologyClaimerOP DensityScoringClaimer::clone() const {
	return new DensityScoringClaimer( *this );
}


void DensityScoringClaimer::add_constraints( core::pose::Pose& ) const {
	// ??
}

void  DensityScoringClaimer::set_defaults() {
	SequenceClaimer::set_defaults();
	anchor_chain_ = ""; //usually anchored to DEFAULT chain
	anchor_residue_ = 0;
}

bool DensityScoringClaimer::read_tag( std::string tag, std::istream& is ) {
	using namespace jumping;
	if ( tag == "anchor" ) {
		is >> anchor_residue_;
		tr.Debug << "Read anchor res = " << anchor_residue_ << std::endl;
	} else if ( SequenceClaimer::read_tag( tag, is ) ) {
		//noop
	} else return false;
	return true;
}

bool DensityScoringClaimer::accept_declined_claim( DofClaim const& was_declined ) {
	tr.Error << "JumpClaim of " << type() << " was declined: " << was_declined << std::endl;
	return false; // full tolerance here ---
}


void DensityScoringClaimer::init_after_reading() {
	if ( !anchor_residue_ ) {
		throw EXCN_Input( "need to specify anchor residue for DensityScoring" );
	}

	set_label( ObjexxFCL::string_of( anchor_residue_ ) + anchor_chain_ );
	if ( !anchor_chain_.size() ) anchor_chain_ = "main"; // ?????

	set_sequence( "X" );
}



} //topology_broker
} //protocols
