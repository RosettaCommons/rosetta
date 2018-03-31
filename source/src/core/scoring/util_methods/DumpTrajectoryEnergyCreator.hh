// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/util_methods/DumpTrajectoryEnergyCreator.hh
/// @brief Creator for an EnergyMethod that dump trajectories to file.
/// @details Dumps trajectories of the minimizer and packer to file when the dump_trajectory
/// ScoreType is enable. Output may be controlled through the dump_trajectory:* flags.
/// @author Brian Coventry (bcov@uw.edu).

#ifndef INCLUDED_core_scoring_util_methods_DumpTrajectoryEnergyCreator_hh
#define INCLUDED_core_scoring_util_methods_DumpTrajectoryEnergyCreator_hh

// Unit header
#include <core/scoring/methods/EnergyMethodCreator.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

// Utility header
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace util_methods {

class DumpTrajectoryEnergyCreator : public core::scoring::methods::EnergyMethodCreator
{
public:

	/// @brief Instantiate a new DumpTrajectoryEnergy.
	///
	virtual core::scoring::methods::EnergyMethodOP create_energy_method( core::scoring::methods::EnergyMethodOptions const &options ) const;

	/// @brief Return the set of score types claimed by the EnergyMethod that
	/// this EnergyMethodCreator creates in its create_energy_method() function.
	virtual ScoreTypes score_types_for_method() const;
};

} // util_methods
} // scoring
} // core

#endif
