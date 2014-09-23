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


#ifndef INCLUDED_protocols_topology_broker_MetalloClaimer_hh
#define INCLUDED_protocols_topology_broker_MetalloClaimer_hh


// Unit Headers
#include <protocols/topology_broker/MetalloClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh>
#include <protocols/topology_broker/FragmentJumpClaimer.hh>
#include <protocols/topology_broker/SequenceClaimer.hh>

#include <protocols/jumping/ResiduePairJumpSetup.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.hh>
//#include <core/fragment/FragSet.hh>
//

// ObjexxFCL Headers

// Utility headers
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
//#include <basic/Tracer.hh>
//#include <basic/options/option.hh>

#include <utility/pointer/ReferenceCount.hh>

//#include <basic/options/option_macros.hh>

//// C++ headers
//#include <fstream>
// AUTO-REMOVED #include <istream>
#include <string>

#include <utility/vector1.hh>


// option key includes


namespace protocols {
namespace topology_broker {

class MetalloClaimer : public SequenceClaimer, public FragmentJumpClaimer {
public:
	MetalloClaimer(); //for factory
	~MetalloClaimer() {};
	//MetalloClaimer( simple_moves::FragmentMoverOP, std::string mover_tag, weights::AbinitioMoverWeightOP weight );
	//MetalloClaimer( simple_moves::FragmentMoverOP );

	virtual TopologyClaimerOP clone() const {
		return TopologyClaimerOP( new MetalloClaimer( *this ) );
	}

	virtual void generate_sequence_claims( claims::DofClaims& dc ) {
		FragmentJumpClaimer::generate_sequence_claims( dc );
		SequenceClaimer::generate_sequence_claims( dc );
	};

	///mainly calls parent function... but is also used to figure out what residue number we are jumping to.
	// virtual void initialize_residues( core::pose::Pose&, claims::SequenceClaimOP init_claim, claims::DofClaims& failed_to_init );

	virtual void generate_claims( protocols::topology_broker::claims::DofClaims& dc);

	///@brief is called after all round1 claims have been approved or retracted -- additional claims can be issued in this round
	//virtual DofClaims finalize_claims( DofClaims& );

 	virtual void initialize_dofs( core::pose::Pose& pose, claims::DofClaims const& init_claims, claims::DofClaims& failed_to_init ) {
 		claims::DofClaims my_failures;
		FragmentJumpClaimer::initialize_dofs( pose, init_claims, my_failures );
 		SequenceClaimer::initialize_dofs( pose, my_failures, failed_to_init );
 	};

	//	virtual bool accept_declined_claim( DofClaim const& was_declined );

	///@brief type() is specifying the output name of the TopologyClaimer
	virtual std::string type() const {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "MetalloClaimer";
	}

	virtual void add_constraints( core::pose::Pose& /*pose*/ ) const;

	//void set_mover( simple_moves::FragmentMoverOP mover ) {
	//		mover_ = mover;
	//	}

// 	void set_mover_tag( std::string const& str ) {
// 		mover_tag_ = str;
// 		if ( mover_ ) mover_->type( str );
// 	}

// 	std::string const& mover_tag() const {
// 		return mover_tag_;
// 	}

//	virtual moves::MoverOP get_mover(	core::pose::Pose const& /*pose*/ ) const;

protected:

	virtual void set_defaults();
	virtual bool read_tag( std::string tag, std::istream& );
	virtual void init_after_reading();

// 	simple_moves::FragmentMover const& mover() const {
// 		return *mover_;
// 	}

//	kinematics::MoveMapOP movemap_;

private:
	jumping::ResiduePairJumpSetupOP jump_setup_;
	jumping::ResiduePairJumpOP residue_pair_jump_;
	std::string ligand_; // if this is ZN the sequence will be Z[ZN]
	core::Size anchor_residue_; //where does this ligand bound to
	std::string anchor_chain_; // a SequenceLabel

	core::Size resolved_anchor_residue_; //residue number of anchor in final pose
}; //class MetalloClaimer

}
}

#endif
