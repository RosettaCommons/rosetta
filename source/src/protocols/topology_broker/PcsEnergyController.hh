// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
/// @file protocols/topology_broker/PcsEnergyController.hh
///
/// @authorv Christophe Schmitz & Oliver Lange
///
////////////////////////////////////////////////

#ifndef INCLUDED_protocols_topology_broker_PcsEnergyController_hh
#define INCLUDED_protocols_topology_broker_PcsEnergyController_hh

// Unit Headers

// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


// ObjexxFCL Headers

// Utility headers

// C++ headers

namespace protocols {
namespace topology_broker {

class PcsEnergyController : public TopologyClaimer {

	typedef TopologyClaimer Parent;

public:

	PcsEnergyController(); // Construct

	~PcsEnergyController() override; // Destruct

	PcsEnergyController(PcsEnergyController const & other); // copy


	TopologyClaimerOP // clone
	clone() const override {
		return TopologyClaimerOP( new PcsEnergyController( *this ) );
	}

	/// @brief type() is specifying the output name of the TopologyClaimer
	std::string
	type() const override {
		return _static_type_name();
	}

	/// @brief return the type name
	static std::string
	_static_type_name() {
		return "PcsEnergyController";
	}

	/// @brief This is called to process each tag
	bool
	read_tag( std::string tag, std::istream & ) override;

	/// @brief Called any time the CLAIMER is being read (before any tag is read)
	void
	set_defaults() override;

	/// @brief This is called each time the stageID is changed
	void
	add_mover(moves::RandomMover& /* random_mover */,
		core::pose::Pose const& /*pose*/,
		abinitio::StageID /*stageID*/, /* abinitio sampler stage */
		core::scoring::ScoreFunction const& /*scorefxn*/,
		core::Real /* progress */ /* progress within stage */
	) override;

	/// @brief Called Each time END CLAIMER is processed
	void
	init_after_reading() override;

private:

};

} // namespace topology_broker
} // namespace protocols

#endif
