// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/UnfoldedStateEnergyCalculator/UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor.fwd.hh
/// @brief  MPI Work Pool Job distributor for UnfoldedStateEnergyCalculator
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_protocols_unfolded_state_energy_calculator_UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor_fwd_hh
#define INCLUDED_protocols_unfolded_state_energy_calculator_UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace unfolded_state_energy_calculator {

class UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor;

typedef utility::pointer::shared_ptr< UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor > UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributorOP;
typedef utility::pointer::shared_ptr< UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor const > UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributorCOP;

} // UnfoldedStateEnergyCalculator
} // protocols

#endif //INCLUDED_protocols_UnfoldedStateEnergyCalculator_UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor_FWD_HH
