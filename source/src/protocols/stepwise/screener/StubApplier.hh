// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/StubApplier.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_StubApplier_HH
#define INCLUDED_protocols_stepwise_screener_StubApplier_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/StubApplier.fwd.hh>
#include <core/kinematics/Stub.hh>

namespace protocols {
namespace stepwise {
namespace screener {

class StubApplier: public StepWiseScreener {

public:

	//constructor
	StubApplier( core::kinematics::Stub & stub );

	//destructor
	~StubApplier();

public:

	virtual
	void
	get_update( sampler::StepWiseSamplerBaseOP sampler );

	virtual
	std::string
	name() const { return "StubApplier"; }

	virtual
	StepWiseScreenerType
	type() const { return STUB_APPLIER; }

private:

	core::kinematics::Stub & stub_;
};

} //screener
} //stepwise
} //protocols

#endif
