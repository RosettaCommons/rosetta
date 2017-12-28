// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/voids_penalty_energy/VoidsPenaltyEnergyCreator.hh
/// @brief Creator for an EnergyMethod intended for packing, which penalizes solutions in which the total volume to fill differs greatly
/// from the total volume of the current set of rotamers.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#ifndef INCLUDED_core_pack_voids_penalty_energy_VoidsPenaltyEnergyCreator_hh
#define INCLUDED_core_pack_voids_penalty_energy_VoidsPenaltyEnergyCreator_hh

// Unit header
#include <core/scoring/methods/EnergyMethodCreator.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

// Utility header
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace voids_penalty_energy {

class VoidsPenaltyEnergyCreator : public core::scoring::methods::EnergyMethodCreator
{
public:

	/// @brief Instantiate a new VoidsPenaltyEnergy.
	///
	virtual core::scoring::methods::EnergyMethodOP create_energy_method( core::scoring::methods::EnergyMethodOptions const &options ) const;

	/// @brief Return the set of score types claimed by the EnergyMethod that
	/// this EnergyMethodCreator creates in its create_energy_method() function.
	virtual core::scoring::ScoreTypes score_types_for_method() const;
};

} // voids_penalty_energy
} // pack
} // core

#endif
