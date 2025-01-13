// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/EvolutionOptions.hh
/// @brief  Class declaration for %EvolutionOptions
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)


#ifndef INCLUDED_protocols_ligand_evolution_EvolutionOptions_HH
#define INCLUDED_protocols_ligand_evolution_EvolutionOptions_HH

// unit headers
#include <protocols/ligand_evolution/EvolutionOptions.fwd.hh>

// package headers
#include <protocols/ligand_evolution/Population.hh>
#include <protocols/ligand_evolution/offspring_factory/OffspringFactory.hh>
#include <protocols/ligand_evolution/Scorer.hh>
#include <protocols/ligand_evolution/selectors/Selector.hh>

// project header
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>

// Utility headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace ligand_evolution {

/// @brief Within EvolutionOptions class combines all options and settings are collected, defaults saved and sanity checks performed before running
class EvolutionOptions {
public:

	/// @brief Prepares the all options with standard settings
	EvolutionOptions();

	/// @brief Loads settings from file and fills missing settings with standards
	explicit EvolutionOptions( const std::string& path_to_option_file );

	// --------------------------
	// Getters
	// --------------------------

	/// @brief Returns maximum number of generations
	core::Size get_max_generations() const;

	/// @brief Returns number of external scoring runs or zero if evolutionary optimization should be used
	core::Size get_external_scoring() const;

	/// @brief Returns all selector names
	utility::vector1< std::string > get_selector_names() const;

	/// @brief Returns the type of one selector
	const std::string& get_selector_type( const std::string& name ) const;

	/// @brief Returns a required parameter for one selector
	core::Real get_selector_parameter( const std::string& name, const std::string& parameter ) const;

	/// @brief Returns all factory names
	utility::vector1< std::string > get_factory_names() const;

	/// @brief Returns the type of one factory
	const std::string& get_factory_type( const std::string& name ) const;

	/// @brief Returns a required parameter for one factory
	core::Real get_factory_parameter( const std::string& name, const std::string& parameter ) const;

	/// @brief Return the list of factory and selector links. The order defines the order of appliance
	const utility::vector1< std::pair< std::string, std::string > >& get_selector_factory_links() const;

	/// @brief Returns the path to a list of smiles which will all be scored
	const std::string& get_path_to_external_smiles() const;

	/// @brief Returns the path to a list of reactions in SMARTS, defining combinatorial rules
	const std::string& get_path_to_reactions() const;

	/// @brief Returns the path to a list of reagents in SMILES defining the fragments for combination
	const std::string& get_path_to_reagents() const;

	/// @brief Returns a map of all used pop initialization procedures
	const std::map< std::string, std::map< std::string, core::Size > >& get_pop_init_options() const;

	/// @brief Returns how many individuals are supported per generation, essentially defining how many can survive
	core::Size get_population_supported_size() const;

	/// @brief Returns the name of the main selector
	const std::string& get_main_selector() const;

	/// @brief Returns the name of the score function used for the main scores within the scorer
	const std::string& get_main_scfx() const;

	/// @brief Returns the path to a score memory file
	const std::string& get_path_score_memory() const;

	/// @brief Returns the name of the ligand chain used by all movers
	const std::string& get_ligand_chain() const;

	/// @brief Returns the number of scoring runs
	core::Size get_n_scoring_runs() const;

	/// @brief Returns the path to directory where all poses will be saved
	const std::string& get_pose_dir_path() const;

	/// @brief Return the similarity penalty for similar ligands in one generation
	core::Real get_similarity_penalty() const;

	/// @brief Return the identity penalty threshold after which the penalty is applied to the lower scoring molecule
	core::Real get_similarity_penalty_threshold() const;

	/// @brief Returns the name of the score term to optimize for
	const std::string& get_main_term() const;

    const std::string& get_protocol_path() const;

    core::pose::PoseOP get_pose_from_stream();

    numeric::xyzVector<core::Real> get_start_xyz() const;

private:

    void parse_cmdline();

	void parse_option_file( const std::string& option_file );

	// TODO write a function which turns all options into string to either print option during start or write default settings to disk

	void check_validity();

	core::Size check_scorer_setup() const;

	core::Size check_selectors() const;

	core::Size check_factories() const;

	core::Size check_pop_init() const;

private:

    /// @brief Counts how many errors are encountered during parsing to allow report of as many errors as possible per program call
    core::Size error_counter_ = 0;

    /// @brief coordinate where centroid of designed ligands will be placed
    numeric::xyzVector<core::Real> xyz_;

    /// @brief Stores all input streams to obtain poses
    core::import_pose::pose_stream::MetaPoseInputStream pose_stream_;

    /// @brief Stores the path to the xml protocol
    std::string protocol_path_;

	/// @brief If set to something else than 0, this will search for smiles list of ligands to score and score them as often as defined here. This is for benchmarking.
	core::Size external_scoring_ = 0;

	/// @brief If external _scoring_ is set to something else than 0, this path will be visited to open a smiles list which will be scored.
	std::string path_to_external_smiles_;

	/// @brief If this is set, the Scorer will try to load scores from this file to use during optimization.
	std::string path_to_score_memory_;

	/// @brief Set the number of generations for evolutionary optimization
	core::Size generations_ = 25;

	/// @brief Path where to find the reaction definitions. The standard settings assume you use a input and run directory structure.
	std::string path_to_reactions_;

	/// @brief Path where to find reagents definitions. The standard settings assume you use a input and run directory structure.
	std::string path_to_reagents_;

	/// @brief Maps selector names to types and parameter list
	std::map< std::string, std::pair< std::string, std::map< std::string, core::Real > > > selector_options_ {
{ "std_elitist", { "elitist", {
// sets how many individuals will be selected
{ "size", 15 },
// sets if selected individuals will be removed from population or remain for further selection
{ "remove", 0 }
} } },
{ "std_tournament", { "tournament", {
// tournament selector selects a subset of individuals to play a "tournament" for being the fittest...
{ "tournament_size", 5 },
// ...and offers them according to rank the victory which they accept with this chance
{ "acceptance_chance", 0.75 },
{ "size", 15 },
{ "remove", 0 }
} } },
// Roulette selector can consider positive values but this distorts chance slightly. 0 maps to false
{ "std_roulette", { "roulette", {
{ "consider_positive", 0 },
{ "size", 15 },
{ "remove", 1 }
} } }
};

	/// @brief Maps offspring factory names to types and parameter list
	std::map< std::string, std::pair< std::string, std::map< std::string, core::Real > > > factory_options_ {
{ "std_mutator", { "mutator", {
{ "size", 30.0 },
{ "reaction_weight", 1.0 },
{ "reagent_weight", 2.0 },
{ "min_similarity", 0.6 },
{ "max_similarity", 0.99 }
} } },
{ "std_crossover", { "crossover", {
{ "size", 30.0 }
} } },
{ "std_identity", { "identity", {
{ "size", 15.0 }
} } }
};

	/// @brief Links a selector to a factory
	utility::vector1< std::pair< std::string, std::string > > selector_factory_links_ {
{ "std_elitist", "std_identity" },
{ "std_roulette", "std_mutator" },
{ "std_roulette", "std_crossover" }
};

	/// @brief Defines how the starting population is generated. For more detail look at %Population
	std::map< std::string, std::map< std::string, core::Size > > pop_init_options_ {
{ "random", {
{ "size", 100 }
} }
};

	/// @brief Defines the size of maximum supported individuals in a population
	core::Size supported_population_size_ = 50;

	/// @brief Sets the name of the main selector. This is used to shrink the population to its supported size.
	std::string main_selector_ = "std_tournament";

	/// @brief Sets how often all movers are applied to a ligand and how often that ligand is scored
	core::Size score_runs_;

	/// @brief Sets the name of the ligand chain
	std::string ligand_chain_;

	/// @brief Sets which score term will be used for optimization
	std::string main_score_term_ = "lid_root2";

	/// @brief If a molecule has better scoring molecules with a Tanimoto similarity exceeding the set threshold, it will receive a score penalty of for each of them. Penalty is added up.
	core::Real similarity_penalty_ = 0.3;

	/// @brief Threshold after which molecules considered being similar and therefore qualified for a penalty.
	core::Real similarity_penalty_threshold_ = 0.95;

	/// @brief Subdirectory within current run directory where all poses will be saved
	std::string pose_dump_directory_;

	/// @brief Name of the scoring function used calculate the final base score
	std::string scoring_function_;

};

}
}

#endif
