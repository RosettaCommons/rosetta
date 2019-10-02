// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/SugarInstantiator.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_SugarInstantiator_HH
#define INCLUDED_protocols_stepwise_screener_SugarInstantiator_HH

#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/screener/SugarInstantiator.fwd.hh>

namespace protocols {
namespace stepwise {
namespace screener {

class SugarInstantiator: public SampleApplier {

public:

	//constructor
	SugarInstantiator( core::pose::Pose & screening_pose,
		core::Size const moving_res,
		core::Distance const o2prime_instantiation_distance_cutoff = 6.0 );

	//destructor
	~SugarInstantiator() override;

public:

	bool
	check_screen() override;

	std::string
	name() const override { return "SugarInstantiator"; }

	StepWiseScreenerType
	type() const override { return SUGAR_INSTANTIATOR; }

	void
	add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover ) override;

private:

	bool check_moving_sugar( core::pose::Pose & pose, core::Size const moving_res );

private:

	core::pose::Pose & screening_pose_;
	core::Size const moving_res_;
	core::Distance const o2prime_instantiation_distance_cutoff_;

	bool instantiate_sugar_;
};

} //screener
} //stepwise
} //protocols

#endif
