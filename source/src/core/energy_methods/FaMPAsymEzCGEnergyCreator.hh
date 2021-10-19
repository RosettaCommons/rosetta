// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/energy_methods/FaMPAsymEzCGCreator.hh
///
/// @brief  Fullatom asymetric EZ potential for CG atoms
/// @details Asymetric EZ potential for CG atoms, from Schramm et al 2012 Structure
///
/// @author  Meghan Franklin (meghanwfranklin@gmail.com)

#ifndef INCLUDED_core_energy_methods_FaMPAsymEzCGEnergyCreator_hh
#define INCLUDED_core_energy_methods_FaMPAsymEzCGEnergyCreator_hh

// Unit Headers
#include <core/scoring/methods/EnergyMethodCreator.hh>

// Package Headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace energy_methods {

/// @brief Asym EZ potential for CG Creator Class
class FaMPAsymEzCGEnergyCreator : public core::scoring::methods::EnergyMethodCreator
{
public:

	/// @brief Instantiate a new FaMPAsymEzCG
	virtual
	core::scoring::methods::EnergyMethodOP
	create_energy_method(
		core::scoring::methods::EnergyMethodOptions const &
	) const;

	/// @brief Return the set of score types claimed by the EnergyMethod
	/// this EnergyMethodCreator creates in its create_energy_method() function
	virtual
	core::scoring::ScoreTypes
	score_types_for_method() const;

};

} // energy_methods
} // core

#endif // INCLUDED_core_energy_methods_FaMPAsymEzCGCreator_hh
