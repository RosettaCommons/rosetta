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


#ifndef INCLUDED_protocols_topology_broker_FibrilModelingClaimer_hh
#define INCLUDED_protocols_topology_broker_FibrilModelingClaimer_hh


// Unit Headers
#include <protocols/topology_broker/FibrilModelingClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <protocols/topology_broker/TopologyClaimer.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <protocols/loops/Loops.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


// option key includes


namespace protocols {
namespace topology_broker {

class FibrilModelingClaimer : public TopologyClaimer {
	typedef TopologyClaimer Parent;

public:

	//c'stor
	FibrilModelingClaimer();
	FibrilModelingClaimer( core::pose::Pose const& input_pose, loops::Loops rigid, int shift = 0 );
	FibrilModelingClaimer( core::pose::Pose const& input_pose, loops::Loops rigid, loops::Loops input_rigid  );

	void
	make_fibril( core::pose::Pose & pose );

	//clone
	TopologyClaimerOP clone() const override {
		return TopologyClaimerOP( new FibrilModelingClaimer( *this ) );
	}

	/// @brief type() is specifying the output name of the TopologyClaimer
	std::string type() const override {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "FibrilModelingClaimer";
	}

	void add_mover(
		moves::RandomMover& random_mover,
		core::pose::Pose const& pose,
		abinitio::StageID stageID,
		core::scoring::ScoreFunction const& scorefxn,
		core::Real progress
	) override;

	void initialize_dofs( core::pose::Pose&,
		claims::DofClaims const& init_claims,
		claims::DofClaims& /*failed_to_init*/ ) override;

	void generate_claims( claims::DofClaims& new_claims ) override;

	/// @brief has to decline foreign BB claims for slave regions
	bool allow_claim( claims::DofClaim const& /*foreign_claim*/ ) override;

protected:

	bool read_tag( std::string tag, std::istream& is ) override;

private:

	/// @brief monomer pose
	core::pose::Pose input_pose_;

	/// @brief regions that can be used for rigid core
	loops::Loops rigid_core_, input_rigid_core_;

	/// @brief align to the monomer before create symmetry
	bool bAlign_;
	int sequence_shift_;

	/// @brief symmetry information
	core::conformation::symmetry::SymmetryInfoOP symminfo_;

};

}
}

#endif
