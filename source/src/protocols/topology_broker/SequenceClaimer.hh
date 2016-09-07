// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_SequenceClaimer_hh
#define INCLUDED_protocols_topology_broker_SequenceClaimer_hh


// Unit Headers
#include <protocols/topology_broker/SequenceClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <protocols/topology_broker/TopologyClaimer.hh>

// Project Headers
#include <core/sequence/Sequence.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
//#include <fstream>
#include <string>

#include <utility/vector1.hh>


// option key includes


namespace protocols {
namespace topology_broker {

class SequenceClaimer : public virtual TopologyClaimer {
	typedef TopologyClaimer Parent;
public:
	SequenceClaimer();
	SequenceClaimer(
		std::string sequence,
		std::string const& label,
		std::string rsd_type_set_identifier
	);

	SequenceClaimerOP shared_from_this() { return utility::pointer::dynamic_pointer_cast<SequenceClaimer>( TopologyClaimer::shared_from_this() ); }

	TopologyClaimerOP clone() const override;

	void generate_sequence_claims( claims::DofClaims& ) override;

	/// @brief is called after all round1 claims have been approved or retracted -- additional claims can be issued in this round
	///if this Sequence has been moved from position 1 --- needs to issue a fixed CUT in the fold-tree
	void generate_claims( claims::DofClaims& ) override;

	// virtual bool allow_claim( DofClaim const& foreign_claim );

	// virtual void initialize_residues( core::pose::Pose&, claims::SequenceClaimOP init_claim, claims::DofClaims& failed_to_init );

	//  virtual void initialize_dofs( core::pose::Pose& pose, claims::DofClaims const& init_claims, claims::DofClaims& failed_to_init );

	/// @brief type() is specifying the output name of the TopologyClaimer
	std::string type() const override {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "SequenceClaimer";
	}

protected:
	bool read_tag( std::string tag, std::istream& ) override;
	void init_after_reading() override;

	void make_sequence_claim();

	void set_sequence( std::string const& str ) {
		input_sequence_ = str;
		sequence_claim_ = nullptr;
	}

	void set_priority( core::Real pr ){
		priority_ =  pr;
	}

private:
	void read_fasta_file( std::string file );

	std::string rsd_type_set_;

	core::Real priority_;
	std::string input_sequence_;


	claims::SequenceClaimOP sequence_claim_;
}; //class SequenceClaimer

}
}

#endif
