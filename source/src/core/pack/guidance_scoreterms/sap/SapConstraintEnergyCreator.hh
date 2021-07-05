// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/sap/SapConstraintEnergyCreator.hh
/// @brief Energy method that enforces the sap_constraint
/// @details sap_constraint constrains your protein to be soluble.
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraintEnergyCreator_hh
#define INCLUDED_core_pack_guidance_scoreterms_sap_SapConstraintEnergyCreator_hh

// Unit header
#include <core/scoring/methods/EnergyMethodCreator.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

// Utility header


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

class SapConstraintEnergyCreator : public core::scoring::methods::EnergyMethodCreator
{
public:

	/// @brief Instantiate a new SapConstraintEnergy.
	///
	core::scoring::methods::EnergyMethodOP create_energy_method( core::scoring::methods::EnergyMethodOptions const &options ) const override;

	/// @brief Return the set of score types claimed by the EnergyMethod that
	/// this EnergyMethodCreator creates in its create_energy_method() function.
	scoring::ScoreTypes score_types_for_method() const override;
};

} //sap
} //guidance_scoreterms
} //pack
} //core

#endif
