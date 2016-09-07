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


#ifndef INCLUDED_protocols_topology_broker_CutBiasClaimer_hh
#define INCLUDED_protocols_topology_broker_CutBiasClaimer_hh


// Unit Headers
#include <protocols/topology_broker/CutBiasClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh>


// Project Headers
//#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
//#include <basic/Tracer.hh>
//#include <basic/options/option.hh>

//#include <basic/options/option_macros.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>

//// C++ headers
//#include <fstream>
#include <string>

#include <utility/vector1.hh>


// option key includes


namespace protocols {
namespace topology_broker {

class CutBiasClaimer : public virtual TopologyClaimer {
public:
	CutBiasClaimer(); //for factory
	CutBiasClaimer( utility::vector1< core::Real > const& );
	CutBiasClaimer( core::fragment::SecondaryStructure const& ss_def );
	CutBiasClaimer( core::fragment::SecondaryStructure const& ss_def, std::string label );

	CutBiasClaimerOP shared_from_this() { return utility::pointer::dynamic_pointer_cast<CutBiasClaimer>( TopologyClaimer::shared_from_this() ); }

	TopologyClaimerOP clone() const override {
		return TopologyClaimerOP( new CutBiasClaimer( *this ) );
	}

	void generate_claims( claims::DofClaims& ) override {};

	/// @brief read definition of Claimer from setup file, i.e., a CLAIMER <type> ... END_CLAIMER block
	void read( std::istream & ) override {};

	/// @brief type() is specifying the output name of the TopologyClaimer
	std::string type() const override {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "CutBiasClaimer";
	}

	void manipulate_cut_bias( utility::vector1< core::Real >& /*cut_bias*/ ) override;

protected:
private:
	utility::vector1< core::Real > cut_bias_;
}; //class CutBiasClaimer

}
}

#endif
