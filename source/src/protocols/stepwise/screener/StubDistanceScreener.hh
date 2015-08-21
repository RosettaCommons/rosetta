// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/StubDistanceScreener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_StubDistanceScreener_HH
#define INCLUDED_protocols_stepwise_screener_StubDistanceScreener_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/StubDistanceScreener.fwd.hh>
#include <core/kinematics/Stub.hh>

namespace protocols {
namespace stepwise {
namespace screener {

class StubDistanceScreener: public StepWiseScreener {

public:

	//constructor
	StubDistanceScreener( core::kinematics::Stub & moving_res_base_stub,
		core::kinematics::Stub const & reference_stub,
		core::Real const max_distance_squared );

	//destructor
	~StubDistanceScreener();

public:

	virtual
	bool check_screen();

	virtual
	std::string
	name() const { return "StubDistanceScreener"; }

	virtual
	StepWiseScreenerType
	type() const { return STUB_DISTANCE; }

	virtual
	void
	fast_forward( sampler::StepWiseSamplerBaseOP sampler );

private:

	core::kinematics::Stub & moving_res_base_stub_;
	core::kinematics::Stub const & reference_stub_;
	core::Real const max_distance_squared_;

};

} //screener
} //stepwise
} //protocols

#endif
