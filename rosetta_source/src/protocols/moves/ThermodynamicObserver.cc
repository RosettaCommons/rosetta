// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/ThermodynamicObserver.cc
/// @brief ThermodynamicObserver methods implemented
/// @author


// Unit Headers
#include <protocols/moves/ThermodynamicObserver.hh>


// Package Headers

// Project Headers
#include <core/pose/Pose.hh>



// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/exit.hh>

// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

//static basic::Tracer TR( "protocols.moves.ThermodynamicObserver" );

namespace protocols {
namespace moves {

///@brief
ThermodynamicObserver::ThermodynamicObserver(
) : Mover()
{
	Mover::type( "ThermodynamicObserver" );
}

ThermodynamicObserver::~ThermodynamicObserver() {}

void
ThermodynamicObserver::initialize_simulation(
	core::pose::Pose & /*pose*/,
	protocols::moves::MetropolisHastingsMover const & /*metropolis_hastings_mover*/
)
{}

void
ThermodynamicObserver::finalize_simulation(
	core::pose::Pose & /*pose*/,
	protocols::moves::MetropolisHastingsMover const & /*metropolis_hastings_mover*/
)
{}

} //moves
} //protocols

