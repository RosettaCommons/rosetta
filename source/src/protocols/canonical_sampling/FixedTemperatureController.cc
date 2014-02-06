// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/canonical_sampling/TemperatureControllerMover.cc
/// @brief TemperatureController methods implemented
/// @author

// Unit Headers
#include <protocols/canonical_sampling/FixedTemperatureController.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Protocols headers
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

namespace protocols {
namespace canonical_sampling {

using namespace std;
using namespace core;
using core::pose::Pose;
using protocols::moves::MoverOP;
using utility::tools::make_vector1;

FixedTemperatureController::FixedTemperatureController(Real temp) 
	: TemperatureController() {

	set_temperatures(make_vector1(temp));
}

FixedTemperatureController::~FixedTemperatureController() {}

MoverOP FixedTemperatureController::clone() const {
	return new FixedTemperatureController(temperature(1));
}

core::Real FixedTemperatureController::temperature_move(
		Pose & pose, MetropolisHastingsMover & mover, Real score) {
	
	return temperature(1);
}

} // canonical_sampling
} // protocols

