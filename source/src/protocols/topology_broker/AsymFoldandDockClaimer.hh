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


#ifndef INCLUDED_protocols_topology_broker_AsymFoldandDockClaimer_hh
#define INCLUDED_protocols_topology_broker_AsymFoldandDockClaimer_hh


// Unit Headers
#include <protocols/topology_broker/AsymFoldandDockClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <protocols/topology_broker/TopologyClaimer.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <protocols/loops/Loops.hh>



namespace protocols {
namespace topology_broker {

class AsymFoldandDockClaimer : public TopologyClaimer {
	typedef TopologyClaimer Parent;

public:

	//c'stor
	AsymFoldandDockClaimer();
	AsymFoldandDockClaimer( core::pose::Pose const& input_pose );

	//clone
	TopologyClaimerOP clone() const override;

	/// @brief type() is specifying the output name of the TopologyClaimer
	std::string type() const override;

	static std::string _static_type_name();

	void add_mover(
		moves::RandomMover& random_mover,
		core::pose::Pose const& pose,
		abinitio::StageID stageID,
		core::scoring::ScoreFunction const& scorefxn,
		core::Real progress
	) override;

	bool read_tag( std::string tag, std::istream& is ) override;

	void initialize_dofs( core::pose::Pose&,
		claims::DofClaims const& init_claims,
		claims::DofClaims& /*failed_to_init*/ ) override;

	void generate_claims( claims::DofClaims& new_claims ) override;

	core::Size docking_jump( core::pose::Pose& pose, core::Size chain_break_res );
private:

	/// @brief starting pose
	core::pose::Pose input_pose_;
	protocols::loops::Loops moving_res_;
	core::Size chain_break_res_, docking_jump_;
	bool docking_local_refine_;

};

}
}

#endif
