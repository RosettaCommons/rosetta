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

#ifndef INCLUDED_protocols_unfolded_state_energy_calculator_UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor_hh
#define INCLUDED_protocols_unfolded_state_energy_calculator_UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor_hh

// Unit headers
#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor.fwd.hh>

// Package headers
#include <protocols/jd2/MPIWorkPoolJobDistributor.hh>
#include <protocols/moves/Mover.fwd.hh>

// Project headers
#include <core/types.hh>

#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <utility/vector1.hh>

#ifdef WIN32
#include <protocols/moves/Mover.hh>
#endif


namespace protocols {
namespace unfolded_state_energy_calculator {

// these function like the mpi tags from MPIWorkPoolJobDistributor.hh
const int UNFOLDED_ENERGY_DATA_TAG( 40 );
const int UNFOLDED_ENERGY_TERMS_TAG( 50 );
const int LENGTH_TLC( 3 );

class UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor : public protocols::jd2::MPIWorkPoolJobDistributor
{
public:
	typedef std::map<std::string, utility::vector1< core::scoring::EMapVector > >::iterator uem_iter;

	/// @brief ctor is protected; singleton pattern
	UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor();

	///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
	///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
	virtual ~UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor();

	/// @brief dummy for master/slave version
	void
	add_unfolded_energy_data( std::string tlc, core::scoring::EMapVector const & scores );

	/// @brief dummy for master/slave version
	void
	set_energy_terms( core::scoring::EMapVector const & weights );

protected:

	/// @brief
	virtual
	void
	master_go( protocols::moves::MoverOP mover );

private:

	/// @brief dummy for master/slave version
	void
	master_add_unfolded_energy_data( std::string tlc, core::scoring::EMapVector const & scores );

	/// @brief dummy for master/slave version
	void
	slave_add_unfolded_energy_data( std::string tlc, core::scoring::EMapVector const & scores );

	/// @brief dummy for master/slave version
	void
	master_set_energy_terms( core::scoring::EMapVector const & weights );

	/// @brief dummy for master/slave version
	void
	slave_set_energy_terms( core::scoring::EMapVector const & weights );

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

#endif //INCLUDED_protocols_UnfoldedStateEnergyCalculator_UnfoldedStateEnergyCalculatorMPIWorkPoolJobDistributor_HH
