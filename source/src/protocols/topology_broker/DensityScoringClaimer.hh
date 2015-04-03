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


#ifndef INCLUDED_protocols_topology_broker_DensityScoringClaimer_hh
#define INCLUDED_protocols_topology_broker_DensityScoringClaimer_hh


// Unit Headers
#include <protocols/topology_broker/DensityScoringClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh>
#include <protocols/topology_broker/SequenceClaimer.hh>
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>


// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.hh>
//#include <core/fragment/FragSet.hh>


// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
//#include <basic/Tracer.hh>
//#include <basic/options/option.hh>

//#include <basic/options/option_macros.hh>

//// C++ headers
//#include <fstream>
#include <string>

#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>


// option key includes


namespace protocols {
namespace topology_broker {

class DensityScoringClaimer : public virtual TopologyClaimer {
public:
	DensityScoringClaimer(); //for factory

	virtual TopologyClaimerOP clone() const;

	///mainly calls parent function... but is also used to figure out what residue number we are jumping to.
	//virtual void initialize_residues( core::pose::Pose&, claims::SequenceClaimOP init_claim, claims::DofClaims& failed_to_init );

	virtual void generate_sequence_claims( claims::DofClaims& );

	virtual void generate_claims( protocols::topology_broker::claims::DofClaims& dc);

	/// @brief type() is specifying the output name of the TopologyClaimer
	virtual std::string type() const {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "DensityScoringClaimer";
	}

	virtual void add_constraints( core::pose::Pose& /*pose*/ ) const;

	virtual bool accept_declined_claim( claims::DofClaim const& was_declined );
protected:

	virtual void set_defaults();
	virtual bool read_tag( std::string tag, std::istream& );
	virtual void init_after_reading();

private:
	core::Size  anchor_residue_;  // where does the VRT jump to
	std::string anchor_chain_;    // a SequenceLabel  << fpd ... is this used??

	core::Size  resolved_anchor_residue_;   // residue number of anchor in final pose
	core::Size  vrt_id_;                    // residue number of VRT in final pose
}; //class MetalloClaimer

}
}

#endif
