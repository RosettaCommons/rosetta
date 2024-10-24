// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/LinearChainbreakEnergyCreator.cc
/// @brief  LinearChainbreakEnergyCreator implementation
/// @author Christopher Miles (cmiles@uw.edu)

#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/energy_methods/LinearChainbreakEnergy.hh>
#include <core/energy_methods/LinearChainbreakEnergyCreator.hh>



namespace core {
namespace energy_methods {


/// @details Creates and initializes a new LinearChainbreakEnergy object with the
/// specified options
core::scoring::methods::EnergyMethodOP
LinearChainbreakEnergyCreator::create_energy_method(const core::scoring::methods::EnergyMethodOptions& opt) const {
	// core::scoring::methods::EnergyMethodOptions is responsible for providing a reasonable default
	// in the event that the caller does not specify a sequence separation
	// constraint. Currently, cst_max_seq_sep() defaults to +inf.
	Size allowable_sequence_separation = opt.cst_max_seq_sep();
	return utility::pointer::make_shared< LinearChainbreakEnergy >(allowable_sequence_separation);
}

core::scoring::ScoreTypes
LinearChainbreakEnergyCreator::score_types_for_method() const {
	core::scoring::ScoreTypes types;
	types.push_back(core::scoring::linear_chainbreak);
	types.push_back(core::scoring::overlap_chainbreak);
	return types;
}

}  // namespace energy_methods
}  // namespace core
