// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/UnfoldedStateEnergyCalculator/UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor.hh
/// @brief  Job distributor for UnfoldedStateEnergyCalculator
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_protocols_unfolded_state_energy_calculator_UnfoldedStateEnergyCalculatorJobDistributor_hh
#define INCLUDED_protocols_unfolded_state_energy_calculator_UnfoldedStateEnergyCalculatorJobDistributor_hh

// Unit headers
#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorJobDistributor.fwd.hh>

// Package headers
#include <protocols/jd2/FileSystemJobDistributor.hh>
#include <protocols/moves/Mover.fwd.hh>

// Project headers
#include <core/types.hh>

#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

// C++ headers
#include <map>
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace unfolded_state_energy_calculator {

class UnfoldedStateEnergyCalculatorJobDistributor : public protocols::jd2::FileSystemJobDistributor
{
public:
	typedef std::map<std::string, utility::vector1< core::scoring::EMapVector > >::iterator uem_iter;

	/// @brief ctor is protected; singleton pattern
	UnfoldedStateEnergyCalculatorJobDistributor();

	///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
	///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
	~UnfoldedStateEnergyCalculatorJobDistributor() override;

	/// @brief
	
	void
	go( protocols::moves::MoverOP mover ) override;

	/// @brief
	void
	add_unfolded_energy_data( std::string tlc, core::scoring::EMapVector const & scores );

	/// @brief
	void
	set_energy_terms( core::scoring::EMapVector const & weights );

private:

	// energy map to hold weights to determin which terms are turned on
	core::scoring::EMapVector energy_terms_;

	// vector of energy maps to hold data
	utility::vector1< core::scoring::EMapVector > unweighted_energies_;

	// map to hold vector of energies for each residue type encountered
	std::map<std::string, utility::vector1< core::scoring::EMapVector > > unweighted_energies_map_;
};

} // UnfoldedStateEnergyCalculator
} // protocols

#endif //INCLUDED_protocols_UnfoldedStateEnergyCalculator_UnfoldedStateEnergyCalculatorJobDistributor_HH
