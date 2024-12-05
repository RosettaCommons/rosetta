// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/EvolutionManager.hh
/// @brief  Class declaration for %EvolutionManager
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)


#ifndef INCLUDED_protocols_ligand_evolution_EvolutionManager_HH
#define INCLUDED_protocols_ligand_evolution_EvolutionManager_HH

// unit headers
#include <protocols/ligand_evolution/EvolutionManager.fwd.hh>
#include <protocols/ligand_evolution/EvolutionOptions.hh>

// package headers
#include <protocols/ligand_evolution/Crossover.hh>
#include <protocols/ligand_evolution/FragmentLibrary.hh>
#include <protocols/ligand_evolution/Mutator.hh>
#include <protocols/ligand_evolution/IdentityFactory.hh>
#include <protocols/ligand_evolution/Population.hh>
#include <protocols/ligand_evolution/Scorer.hh>
#include <protocols/ligand_evolution/WorkManager.hh>


namespace protocols {
namespace ligand_evolution {

/// @brief The EvolutionManager combines all required resources for the evolutionary ligand optimization and handles them as needed.
class EvolutionManager {
public:

	/// @brief If rank is different than 0 a mpi system is expected
	explicit EvolutionManager( int rank );

	/// @brief Main function that handles evolutionary optimization
	void run( int mpi_size = 1 );

	void init();

	/// @brief Prints scores of all individuals
	std::string print_scores() const;

private:

	/// @brief Scores the entire population
	void score();

	void init_mpi();

	/// @brief Calculates the current quantiles
	void calculate_quantiles();

	/// @brief Saves information about the population and its history to disk. For debugging and benchmarking purposes
	void write_population_information() const;

private:

	/// @brief Defines after how many generations optimization will stop
	core::Size max_generations_;

	/// @brief Contains and manages all evolutionary information for all ligands
	Population population_;

	/// @brief Manages the scoring of ligands as well as keeping track of scores
	ScorerOP scorer_ = nullptr;

	/// @brief Handles mpi communication to distribute the scoring work
	WorkManagerOP work_manager_ = nullptr;

	/// @brief Holds all molecular information about the protein and the ligands
	FragmentLibrary library_;

	/// @brief The selector which will be applied to advance the generation
	core::Size main_selector_;

	/// @brief Holds all selectors
	utility::vector1< SelectorOP > selectors_;

	/// @brief Maps selector names to indices in selector list
	std::map< std::string, core::Size > selector_map_;

	/// @brief Holds all offspring factories
	utility::vector1< OffspringFactoryOP > factories_;

	/// @brief Maps factory names to indices in factory list
	std::map< std::string, core::Size > factory_map_;

	/// @brief Holds arrays of size 5, defining index of selector, index of factory, selection size, offspring size, if parent will be kept in pool
	/// @details Factories will be called in order of this array
	utility::vector1< utility::vector1< core::Size > > offspring_options_;

	int rank_ = 0;

	/// @brief Stores statistical information about score distribution within the population
	utility::vector1< core::Real > quantiles = { 9999.9, 9999.9, 9999.9, 9999.9 };

	/// @brief Used to distinguish between scoring of external molecules and internal evolutionary optimization
	core::Size external_scoring_;

};

}
}

#endif
