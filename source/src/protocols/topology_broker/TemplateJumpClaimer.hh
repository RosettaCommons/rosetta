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


#ifndef INCLUDED_protocols_topology_broker_TemplateJumpClaimer_hh
#define INCLUDED_protocols_topology_broker_TemplateJumpClaimer_hh


// Unit Headers
#include <protocols/topology_broker/TemplateJumpClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/FragmentJumpClaimer.hh>
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <protocols/topology_broker/weights/AbinitioMoverWeight.hh>
#include <protocols/abinitio/Templates.hh>
// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/JumpSample.hh>
#include <core/scoring/dssp/PairingsList.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <protocols/jumping/SheetBuilder.hh>
// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
//#include <basic/Tracer.hh>
//#include <basic/options/option.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


//#include <basic/options/option_macros.hh>

//// C++ headers
//#include <fstream>


// option key includes


namespace protocols {
namespace topology_broker {

/// @brief hacky wrapper to keep the old Template code alive a bit longer
/// this claimer deals with the Jumpy part of the Templates.
class TemplateJumpClaimer : public FragmentJumpClaimer {
public:
	TemplateJumpClaimer(); //for factory
	TemplateJumpClaimer( std::string config_file,  weights::AbinitioMoverWeightOP weight = NULL );

	virtual TopologyClaimerOP clone() const {
		return TopologyClaimerOP( new TemplateJumpClaimer( *this ) );
	}

	void read_config_file( std::string const& file );
	void read_topol_file( std::string const& file );

	/// @brief type() is specifying the output name of the TopologyClaimer
	virtual std::string type() const {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "TemplateJumpClaimer";
	}

protected:
	virtual bool read_tag( std::string tag, std::istream& is );
	virtual void init_after_reading();
private:
	// info about homologues structures --- if available
	abinitio::TemplatesOP templates_;

	// or should we make jumps in the old-fashioned way with SheetBuilder?
	core::scoring::dssp::PairingsList pairings_;
	bool bRandomSheet_;
	core::fragment::SecondaryStructureOP ss_def_;
	jumping::SheetBuilder::SheetTopology sheets_;

}; //class TemplateJumpClaimer

}
}

#endif
