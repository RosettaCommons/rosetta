// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/UnfoldedStateEnergyCalculator/UnfoldedStateEnergyCalculatorMover.fwd.hh
/// @brief forward declaration of UnfoldedStateEnergyCalculatorMover
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_protocols_unfolded_state_energy_calculator_UnfoldedStateEnergyCalculatorMover_fwd_hh
#define INCLUDED_protocols_unfolded_state_energy_calculator_UnfoldedStateEnergyCalculatorMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace unfolded_state_energy_calculator {

class UnfoldedStateEnergyCalculatorMover;

typedef utility::pointer::shared_ptr< UnfoldedStateEnergyCalculatorMover > UnfoldedStateEnergyCalculatorMoverOP;
typedef utility::pointer::shared_ptr< UnfoldedStateEnergyCalculatorMover const > UnfoldedStateEnergyCalculatorMoverCOP;

} // UnfoldedStateEnergyCalculator
} // protocols

#endif //INCLUDED_protocols_UnfoldedStateEnergyCalculator_UnfoldedStateEnergyCalculatorMover_FWD_HH
