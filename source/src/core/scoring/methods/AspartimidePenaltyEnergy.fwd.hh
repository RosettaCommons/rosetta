// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/AspartimidePenaltyEnergy.fwd.hh
/// @brief  Forward declarations for the AspartimidePenaltyEnergy.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#ifndef INCLUDED_core_scoring_methods_AspartimidePenaltyEnergy_fwd_hh
#define INCLUDED_core_scoring_methods_AspartimidePenaltyEnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class AspartimidePenaltyEnergy;
typedef  utility::pointer::weak_ptr< AspartimidePenaltyEnergy > AspartimidePenaltyEnergyAP;
typedef  utility::pointer::weak_ptr< AspartimidePenaltyEnergy const > AspartimidePenaltyEnergyCAP;
typedef  utility::pointer::shared_ptr< AspartimidePenaltyEnergy > AspartimidePenaltyEnergyOP;
typedef  utility::pointer::shared_ptr< AspartimidePenaltyEnergy const > AspartimidePenaltyEnergyCOP;

} // namespace methods
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_methods_AspartimidePenaltyEnergy_HH
