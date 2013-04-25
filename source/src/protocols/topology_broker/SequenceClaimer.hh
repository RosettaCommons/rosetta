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


#ifndef INCLUDED_protocols_topology_broker_SequenceClaimer_hh
#define INCLUDED_protocols_topology_broker_SequenceClaimer_hh


// Unit Headers
#include <protocols/topology_broker/SequenceClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/DofClaim.fwd.hh>
#include <protocols/topology_broker/TopologyClaimer.hh>

// Project Headers

// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
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

class SequenceClaimer : public virtual TopologyClaimer {
	typedef TopologyClaimer Parent;
public:
	SequenceClaimer();
	SequenceClaimer( std::string const& sequence, std::string const& rsd_type_set_identifier, std::string label );

	virtual TopologyClaimerOP clone() const;

	virtual void generate_sequence_claims( DofClaims& );

	///@brief is called after all round1 claims have been approved or retracted -- additional claims can be issued in this round
	///if this Sequence has been moved from position 1 --- needs to issue a fixed CUT in the fold-tree
	virtual void generate_claims( DofClaims& );

	//	virtual bool allow_claim( DofClaim const& foreign_claim );

	virtual void initialize_residues( core::pose::Pose&, SequenceClaimOP init_claim, DofClaims& failed_to_init );

 	virtual void initialize_dofs( core::pose::Pose& pose, DofClaims const& init_claims, DofClaims& failed_to_init );

// 	virtual bool reinitialize_residues() {
// 		return offset_ == 0;
// 	}
	///@brief type() is specifying the output name of the TopologyClaimer
	virtual std::string type() const {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "SequenceClaimer";
	}

	core::Size offset() {
		return offset_; //where in the pose does this sequence start?
	}

	void set_sequence( std::string const& str );

protected:
	virtual bool read_tag( std::string tag, std::istream& );
	virtual void init_after_reading();

private:
	std::string sequence_;
	std::string annotated_sequence_;
	std::string rsd_type_set_;
	core::Size offset_;
	core::Size nr_res_;
}; //class SequenceClaimer

}
}

#endif
