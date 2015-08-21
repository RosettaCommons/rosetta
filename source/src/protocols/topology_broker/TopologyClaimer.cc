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
#include <protocols/topology_broker/TopologyClaimer.hh>

// Package Headers
#include <protocols/topology_broker/Exceptions.hh>

// Project Headers
#include <protocols/topology_broker/weights/LargeFragWeight.hh>
#include <protocols/topology_broker/weights/SmoothFragWeight.hh>
#include <protocols/topology_broker/weights/SmallFragWeight.hh>
#include <protocols/topology_broker/weights/ConstAbinitioMoverWeight.hh>
#include <protocols/moves/MoverContainer.hh>


// Utility headers
#include <basic/Tracer.hh>


// C/C++ headers
#ifdef WIN32
#include <algorithm>
#include <iterator>
#endif

#include <vector>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

/// @details Auto-generated virtual destructor
//TopologyClaimer::~TopologyClaimer();

using namespace core;

void TopologyClaimer::initialize_dofs( core::pose::Pose&, claims::DofClaims const& init_dofs, claims::DofClaims& failed_to_init ) {
	claims::DofClaims my_claims;
	for ( claims::DofClaims::const_iterator it = init_dofs.begin(), eit = init_dofs.end();
			it != eit; ++it ) {
		if ( (*it)->owner().lock().get() == this ) {
			my_claims.push_back( *it );
		}
	}
	if ( my_claims.size() ) {
		tr.Warning << "[WARNING]" << type() << "did not initialize dofs as requested" << std::endl;
	}
	std::copy( my_claims.begin(), my_claims.end(), std::back_inserter( failed_to_init ) );
}

/*void TopologyClaimer::initialize_residues( core::pose::Pose&, claims::SequenceClaimOP init_claim, claims::DofClaims& failed_to_init ) {
runtime_assert( init_claim->owner()==this );
failed_to_init.push_back( init_claim );
tr.Warning << "[WARNING]" << type() << "did not initialize residues as requested for claim..." << *init_claim << std::endl;
}*/

void TopologyClaimer::add_mover(
	moves::RandomMover& random_mover,
	core::pose::Pose const& pose,
	abinitio::StageID stageID, /* abinitio sampler stage */
	core::scoring::ScoreFunction const& /*scorefxn*/, /* scorefxn of this stage */
	core::Real progress /* progress within stage */
) {
	moves::MoverOP mover = get_mover( pose );
	if ( mover ) {
		//if(tr.Debug.visible()){tr.Debug << "current mover is:  " << mover->type() << std::endl;}
		random_mover.add_mover( mover, abinitio_mover_weight_->weight( stageID, progress ) );
	}
}

void TopologyClaimer::read_mover_weight( std::istream& is ) {
	std::string type;
	core::Real weight;
	is >> type >> weight;
	if ( type == "LargeStage" ) { //sample when ClassicAbinitio likes to have "large" fragments
		set_mover_weight( weights::AbinitioMoverWeightOP( new weights::LargeFragWeight( weight ) ) );
	} else if ( type == "SmallStage" ) { //sample when ClassicAbinitio likes to have "small" fragments
		set_mover_weight( weights::AbinitioMoverWeightOP( new weights::SmallFragWeight( weight ) ) );
	} else if ( type == "SmoothStage" ) { //sample when ClassicAbinitio likes to have "small" fragments
		set_mover_weight( weights::AbinitioMoverWeightOP( new weights::SmoothFragWeight( weight ) ) );
	} else if ( type == "AllStage" ) { //always on
		set_mover_weight( weights::AbinitioMoverWeightOP( new weights::ConstAbinitioMoverWeight( weight ) ) );
	} else {
		throw EXCN_Input( "weight can only by one of LargeStage, SmallStage or AllStage " );
	}
}

void TopologyClaimer::unknown_tag( std::string tag, std::istream& is ) const {
	std::string line; getline(is, line);
	throw EXCN_Input ("ERROR reading broker-setup. unknown tag: " + tag + " in line:\n " + tag + line + "\n");
}

void TopologyClaimer::read( std::istream& is ) {
	set_defaults();
	std::string tag;
	while ( is >> tag && tag != "END_CLAIMER" ) {
		if ( tag[ 0 ] == '#' ) {
			getline( is, tag );
			continue;
		}
		tr.Trace << "READ_SETUP (" << type() << "): tag = " << tag << std::endl;
		if ( !read_tag( tag, is ) ) unknown_tag( tag, is );
	}
	init_after_reading();
}

bool TopologyClaimer::read_tag( std::string tag, std::istream& is ) {
	if ( tag == "LABEL" ) {
		is >> label_;
		return true;
	}

	std::string arg;
	is >> arg;
	tr.Error << "[ERROR]: The tag '" << tag << "' with argument '" << arg << "' was not recognized." << std::endl;
	return false;
}


void TopologyClaimer::set_defaults() {
	label_ = "DEFAULT";
}


} //topology_broker
} //protocols
