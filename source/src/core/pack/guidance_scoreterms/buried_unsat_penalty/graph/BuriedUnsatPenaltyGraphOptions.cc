// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphOptions.cc
/// @brief Options container for the BuriedUnsatPenaltyGraph.  Initialized by the BuriedUnsatPenalty energy method.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphOptions.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "core.pack.guidance_scoreterms.buried_unsat_penalty.graph.BuriedUnsatPenaltyGraphOptions" );


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace buried_unsat_penalty {
namespace graph {

/// @brief Options constructor
BuriedUnsatPenaltyGraphOptions::BuriedUnsatPenaltyGraphOptions(
	core::Real const angle_exponent,
	core::Real const angle_shift_factor,
	core::Real const dist_exponent,
	core::Real const dist_midpoint,
	core::Real const burial_threshold,
	core::Real const hbond_energy_threshold
):
	utility::pointer::ReferenceCount(),
	angle_exponent_(angle_exponent),
	angle_shift_factor_(angle_shift_factor),
	dist_exponent_(dist_exponent),
	dist_midpoint_(dist_midpoint),
	burial_threshold_(burial_threshold),
	hbond_energy_threshold_(hbond_energy_threshold)
{}

/// @brief Default destructor
BuriedUnsatPenaltyGraphOptions::~BuriedUnsatPenaltyGraphOptions(){}

/// @brief Clone method: return a copy of the original object, by owning pointer.
BuriedUnsatPenaltyGraphOptionsOP
BuriedUnsatPenaltyGraphOptions::clone() const {
	return BuriedUnsatPenaltyGraphOptionsOP( new BuriedUnsatPenaltyGraphOptions( *this ) );
}

/// @brief Print information about this object.
void
BuriedUnsatPenaltyGraphOptions::show(
	std::ostream & out
) const {
	out << "angle_exponent = " << angle_exponent_ << std::endl;
	out << "angle_shift_factor = " << angle_shift_factor_ << std::endl;
	out << "dist_exponent = " << dist_exponent_ << std::endl;
	out << "dist_midpoint = " << dist_midpoint_ << std::endl;
	out << "burial_threshold = " << burial_threshold_ << std::endl;
	out << "hbond_energy_threshold = " << hbond_energy_threshold_ << std::endl;
}

} //core
} //pack
} //guidance_scoreterms
} //buried_unsat_penalty
} //graph






