// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/UnfoldedStateEnergyCalculator/UnfoldedStateEnergyCalculatorUtil.hh
/// @brief  Utility functions common to both UnfoldedStateEnergyCalculatorjob distributors
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_protocols_unfolded_state_energy_calculator_UnfoldedStateEnergyCalculatorUtil_hh
#define INCLUDED_protocols_unfolded_state_energy_calculator_UnfoldedStateEnergyCalculatorUtil_hh

// Project headers
#include <core/types.hh>

#include <core/scoring/EnergyMap.fwd.hh>

// Utility headers
#include <utility/vector1.hh>




namespace protocols {
namespace unfolded_state_energy_calculator {

///@brief
void
calc_all_averages( utility::vector1< core::scoring::EMapVector > unweighted_energies, core::scoring::EMapVector energy_terms );

///@brief
core::Real
calc_vector_mean( utility::vector1< core::Real> const & data );

///@brief
core::Real
calc_vector_median( utility::vector1< core::Real> const & data );

///@brief
core::Real
calc_vector_mode( utility::vector1< core::Real> const & data );

///@brief
core::Real
calc_vector_boltzmann( utility::vector1< core::Real> const & data );

} // UnfoldedStateEnergyCalculator
} // protocols

#endif //INCLUDED_protocols_UnfoldedStateEnergyCalculator_UnfoldedStateEnergyCalculatorUtil_HH
