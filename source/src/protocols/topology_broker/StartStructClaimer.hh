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


#ifndef INCLUDED_protocols_topology_broker_StartStructClaimer_hh
#define INCLUDED_protocols_topology_broker_StartStructClaimer_hh


// Unit Headers
#include <protocols/topology_broker/StartStructClaimer.fwd.hh>
#include <protocols/topology_broker/FragmentClaimer.hh>

// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/loops/Loops.hh>


//#include <core/fragment/FragSet.hh>


// ObjexxFCL Headers

// Utility headers


//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
//#include <basic/Tracer.hh>
//#include <basic/options/option.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <set>

//#include <basic/options/option_macros.hh>

//// C++ headers
//#include <fstream>
#include <string>

#include <utility/vector1.hh>


// option key includes


namespace protocols {
namespace topology_broker {

class StartStructClaimer : public virtual FragmentClaimer {
public:
	StartStructClaimer(); //for factory
	StartStructClaimer( core::pose::Pose const& /*idealized*/ );

	StartStructClaimerOP shared_from_this() { return utility::pointer::dynamic_pointer_cast<StartStructClaimer>( FragmentClaimer::shared_from_this() ); }

	TopologyClaimerOP clone() const override {
		return TopologyClaimerOP( new StartStructClaimer( *this ) );
	}

	/// @brief type() is specifying the output name of the TopologyClaimer
	std::string type() const override {
		return _static_type_name();
	}

	/// @brief overloaded to allow perturbation of start structure
	void initialize_dofs( core::pose::Pose&, claims::DofClaims const& init_claims, claims::DofClaims& failed_to_init ) override;
	void generate_claims( claims::DofClaims& ) override;


	static std::string _static_type_name() {
		return "StartStructClaimer";
	}

	void new_decoy( core::pose::Pose const& ) override;
	void new_decoy( ) override;
	bool read_tag( std::string tag, std::istream & ) override;

	moves::MoverOP get_mover( core::pose::Pose const& /*pose*/ ) const override
	{ return nullptr; }; /*does not provide mover*/

protected:
	void generate_init_frags( core::pose::Pose const& );
private:

	/// @brief use the job input pose to get starting structure
	bool bUseInputPose_;

	/// @brief perturb start torsions by gaussian()*perturb_
	core::Real perturb_;

	core::pose::Pose start_pose_;
}; //class StartStructClaimer


}
}

#endif
