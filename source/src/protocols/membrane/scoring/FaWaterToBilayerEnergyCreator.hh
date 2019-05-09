// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/membrane/scoring/FaWaterToBilayerEnergyCreator.hh
/// @brief Water to Bilayer transfer free energy for atom groups derived from
/// Moon & Fleming thermodynamic side chain hydrophobicity scale
/// @author  Rebecca Alford (ralford3@jhu.edu)

#ifndef INCLUDED_protocols_membrane_scoring_FaWaterToBilayerEnergyCreator_hh
#define INCLUDED_protocols_membrane_scoring_FaWaterToBilayerEnergyCreator_hh

// Unit Headers
#include <core/scoring/methods/EnergyMethodCreator.hh>

// Package Headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace membrane {
namespace scoring {

/// @brief Creator Class for Membrane CBeta Energy Method
class FaWaterToBilayerEnergyCreator : public core::scoring::methods::EnergyMethodCreator {

public:

	/// @brief Instantiate a new MPEnvEnergy
	virtual
	core::scoring::methods::EnergyMethodOP
	create_energy_method(
		core::scoring::methods::EnergyMethodOptions const &
	) const;

	/// @brief Return MPEnv Score Type Claimed by this energy method
	virtual
	core::scoring::ScoreTypes
	score_types_for_method() const;
};


} // scoring
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_scoring_FaWaterToBilayerEnergyCreator_hh
