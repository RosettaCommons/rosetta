// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/canonical_sampling/ThermodynamicMover.cc
/// @brief ThermodynamicMover methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/ThermodynamicMover.hh>


// Package Headers

// Project Headers


// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/exit.hh>

#include <core/id/DOF_ID_Range.hh>
#include <utility/vector1.hh>


// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

//static basic::Tracer TR( "protocols.canonical_sampling.ThermodynamicMover" );

namespace protocols {
namespace canonical_sampling {

/// @brief
ThermodynamicMover::ThermodynamicMover(
) : Mover()
{
	Mover::type( "ThermodynamicMover" );
}

ThermodynamicMover::~ThermodynamicMover() {}

void
ThermodynamicMover::initialize_simulation(
	core::pose::Pose & /*pose*/,
	protocols::canonical_sampling::MetropolisHastingsMover const &, /*metropolis_hastings_mover*/
	core::Size //default=0; non-zero if trajectory is restarted
)
{}

core::Real
ThermodynamicMover::last_proposal_density_ratio()
{
	return 1;
}

void
ThermodynamicMover::observe_after_metropolis(
	protocols::canonical_sampling::MetropolisHastingsMover const & /*metropolis_hastings_mover*/
)
{}

void
ThermodynamicMover::finalize_simulation(
	core::pose::Pose & /*pose*/,
	protocols::canonical_sampling::MetropolisHastingsMover const & /*metropolis_hastings_mover*/
)
{}

bool
ThermodynamicMover::is_multi_trial()
{
	return false;
}

core::Real
ThermodynamicMover::last_inner_score_delta_over_temperature()
{
	return 0;
}

protocols::canonical_sampling::MetropolisHastingsMoverAP
ThermodynamicMover::metropolis_hastings_mover()
{
	return protocols::canonical_sampling::MetropolisHastingsMoverAP();
}

void
ThermodynamicMover::set_metropolis_hastings_mover(
	protocols::canonical_sampling::MetropolisHastingsMoverAP //metropolis_hastings_mover
)
{}


/// @details Torsion DOFs that would be returned by torsion_id_ranges() are not
/// returned by this function too.
utility::vector1<core::id::DOF_ID_Range>
ThermodynamicMover::dof_id_ranges(
	core::pose::Pose & //pose
)
{
	return utility::vector1<core::id::DOF_ID_Range>();
}

} //moves
} //protocols

