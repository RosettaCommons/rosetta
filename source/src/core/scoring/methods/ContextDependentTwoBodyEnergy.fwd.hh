// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh
/// @brief  Forward declaration for the Context-Dependent (ShortRange) Two-Body EnergyMethod base class.
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_methods_ContextDependentTwoBodyEnergy_fwd_hh
#define INCLUDED_core_scoring_methods_ContextDependentTwoBodyEnergy_fwd_hh

// Unit headers

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class ContextDependentTwoBodyEnergy;
typedef utility::pointer::shared_ptr< ContextDependentTwoBodyEnergy > ContextDependentTwoBodyEnergyOP;
typedef utility::pointer::shared_ptr< ContextDependentTwoBodyEnergy const > ContextDependentTwoBodyEnergyCOP;

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_ScoreFunction_HH
