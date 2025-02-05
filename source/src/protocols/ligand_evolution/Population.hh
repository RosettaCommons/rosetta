// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/Population.hh
/// @brief  Declaration of the %Population class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

#ifndef INCLUDED_protocols_ligand_evolution_Population_HH
#define INCLUDED_protocols_ligand_evolution_Population_HH

// unit headers
#include <protocols/ligand_evolution/Population.fwd.hh>

// package headers
#include <protocols/ligand_evolution/Individual.hh>
#include <protocols/ligand_evolution/selectors/Selector.hh>
#include <protocols/ligand_evolution/EvolutionOptions.hh>
#include <protocols/ligand_evolution/Scorer.hh>


// utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace ligand_evolution {

/// @ Summarizes and contains all Individuals of one generation
class Population {
public:

	Population() = default;
	~Population() = default;

	/// @brief initializes setting from options object
	void initialize_from_evotoptions( EvolutionOptions const& options, FragmentLibrary const& library, Scorer const& scorer );

	/// @brief Exposes the individuals vector for const iterator access
	utility::vector1< Individual > const& individuals() const;

	/// @brief Expose the individuals vector for normal iterator access
	utility::vector1< Individual >& individuals();

	/// @brief Sorts all individuals depending on their score
	void sort();

	/// @brief Returns true if the individuals are sorted with lowest (and therefore best) score first
	bool is_sorted();

	/// @brief Returns the size of this population
	core::Size size() const;

	/// @brief Sets the supported size of this population. It will get reduced to this size through selective pressure
	void set_supported_size( core::Size supported_size );

	/// @brief Access specific individual
	Individual const& individual( core::Size index ) const;

	/// @brief Non-const access specific individual
	Individual & individual( core::Size index );

	/// @brief Adds random individuals to a population
	void add_random( core::Size n_random_individuals, FragmentLibrary const& lib );

	/// @brief Adds new (and potentially unscored individuals to the population)
	void add_individuals( utility::vector1< Individual > const& new_individuals );

	/// @brief Adds new individuals based on provided LigandIdentifiers
	void add_individuals( utility::vector1< LigandIdentifier > const& new_individuals );

	/// @brief Removes unselected individuals, raises the generation counter and sorts individuals
	void next_generation( Selector const& selector );

	/// @brief Adds an individual to this population and sets its id
	void add_individual( Individual const& individual );

	/// @brief Returns the current generation
	core::Size generation() const;

	/// @brief Removes individuals with the given index and returns them
	utility::vector1< Individual > remove_individuals( utility::vector1< core::Size > indices );

	/// @brief Returns a copy of selected individuals
	utility::vector1< Individual > individuals( utility::vector1< core::Size > const& indices );

	/// @brief Replaces the internal population with a new one
	void replace_population( utility::vector1< Individual > const& individuals );

	/// @brief For debugging and benchmarking purpose to observe population development
	utility::vector1< utility::vector1< core::Size > > const& expose_generation_log() const;

	/// @brief For debugging and benchmarking purpose to observe population development
	utility::vector1< std::pair< core::Size, core::Size > > const& expose_inheritance_graph() const;

private:

	/// @brief Updates and returns the sorting state. Shouldn't be needed to be called usually.
	bool check_sorting();

	/// @brief Saves ids of all sorted individuals
	void log_generation();

	/// @brief Saves the data of the given individual in the inheritance graph
	void expand_inheritance_graph( Individual const& individual );

private:

	utility::vector1< Individual > individuals_;

	core::Size generation_ = 0;

	bool sorting_guaranteed_ = false;

	/// @brief Keeps track of the next available and unused id
	core::Size next_id_ = 1;

	core::Size supported_size_ = 0;

	utility::vector1< utility::vector1< core::Size > > generation_log_;

	/// @brief For debugging and benchmarking purposes. Each entry represents a node in an inheritance graph
	///        First is child, second is the parent. If more than one parent is present, more than one entry is created.
	utility::vector1< std::pair< core::Size, core::Size > > inheritance_graph_;

};


}
}

#endif
