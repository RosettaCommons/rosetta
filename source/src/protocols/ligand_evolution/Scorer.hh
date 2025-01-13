// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/Scorer.hh
/// @brief  Class declaration of the %Scorer class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

#ifndef INCLUDED_protocols_ligand_evolution_Scorer_HH
#define INCLUDED_protocols_ligand_evolution_Scorer_HH

// unit headers
#include <protocols/ligand_evolution/Scorer.fwd.hh>

// project headers
#include <protocols/ligand_evolution/FragmentLibrary.hh>
#include <protocols/ligand_evolution/Population.hh>
#include <protocols/moves/Mover.hh>

// utility headers
#include <utility/vector1.hh>

// C/C++ headers
#include <set>

namespace protocols {
namespace ligand_evolution {

/// @brief: The %Scorer computes the score of each %Individual in a %Population. It applies a list of movers, calculates a variety of score terms and combines them.
/// The scorer also collects statistics and summarizes information.
class Scorer {
public:

	Scorer( FragmentLibrary& library, core::Size n_runs, char ligand_chain );

	/// @brief Sets the score for all individuals in the population, calculates it where needed
	/// @details If the same ligand is multiple times present within the population, each occures after the first gets an increasing penalty.
	void score_population( Population& pop );

	/// @brief Checks if this individuals ligand was already scores, scores it if not, and sets the score for the individual
	void score_individual( Individual& individual );

	/// @brief Sets the current ligand and calls all steps to score it completely
	void score_ligand( LigandIdentifier const& ligand, std::string const& smiles = "", core::Size save_n_scores = 0 );

	/// @brief Adds a mover
	void add_mover( moves::MoverOP const& mover );

	/// @brief Sets the current ligand on which should be worked
	void set_ligand( LigandIdentifier const& ligand, std::string const& smiles = "" );

	/// @brief Checks if a ligand is currently set
	bool has_ligand() const;

	/// @brief Performs the next scoring step on the current individual. Returns true, if this was the last step
	bool next_step( core::Size save_n_scores = 0 );

	/// @brief Returns true if the ligand was already scored
	bool is_scored( LigandIdentifier const& ligand ) const;

	/// @brief Returns all scores in a map with a string identifier
	std::map< std::string, core::Real > const& get_scores( LigandIdentifier const& ligand ) const;

	/// @brief Sets which term should be used for optimization
	void set_main_term( std::string const& score_term );

	/// @brief Sets the penalty for ligands which are appear more than once within a population
	/// @details The penalty is a relative fraction of the current main score for each ligand. It starts at 0 and increases with each additional ligand by the base penalty
	void set_base_similarity_penalty( core::Real base_penalty );

	/// @brief Sets the path where poses should be saved. Leave it empty to save in the current working directory
	void set_pose_path( std::string const& path );

	/// @brief Sets the desired score function
	void set_score_function( core::scoring::ScoreFunctionOP score_function );

	/// @brief Sets all score terms for an identifier if not set yet.
	/// @details This functions expects that the scores were calculated by another instance of a Scorer which took care of all runs and wrote the best pose to disk
	void set_scores( LigandIdentifier const& identifier, std::map< std::string, core::Real > const& scores );

	/// @brief Return how many score terms are calculated
	core::Size n_score_terms() const;

	/// @brief Returns the raw scores of an identifier for mpi communication
	utility::vector1< core::Real > get_raw_scores( LigandIdentifier const& identifier ) const;

	/// @brief Loads scores for combinations of reagents and reactions from file. These are used during the run instead of rescoring. Format is the output format of RLE.
	void load_scores( std::string const& path );

	/// @brief Transforms raw double scores into a proper map and calls set_scores
	/// @detail Expects scores in the same order as score_terms
	void set_raw_scores( LigandIdentifier const& identifier, double const* raw_scores );

	/// @brief Writes all results to disk
	void save_results() const;

	/// @brief Writes multiple scores per ligand to disk
	void save_external_scoring_results( core::Size rank );

	/// @brief For debugging and benchmarking purposes
	std::map< LigandIdentifier, std::set< core::Size > > const& expose_id_memory() const;

	/// @brief Checks if score for this ligand is available in memory and loads it if possible.
	bool check_memory( LigandIdentifier const& ligand );

	/// @brief Returns a list of LigandIdentifiers for loaded scores sorted by their main term score
	/// @param size Sets how many identifiers should be returned. If <= 0, all will be returned
	utility::vector1< LigandIdentifier > get_best_loaded( core::Size size ) const;

	void set_similarity_penalty_threshold( core::Real threshold );

private:

	/// @brief Creates a pose and rotamers for the current ligand
	void create_pose();

	/// @brief Applies all movers in order to the current ligands pose
	void apply_movers();

	/// @brief Calculates all scores for the current ligand and saves them if lower than the old ones
	/// @details If save_all_scores is set, all scores are save in memory. This should only be used for external smiles scoring since it overrides a lot of otherwise needed
	///          and/or desired features.
	void calculate_scores( core::Size save_n_scores = 0 );

	/// @brief Dumps the current best pose to disk with the ligand identifier as name
	void dump_pose();

	utility::vector1< std::string > ligand_line( LigandIdentifier const& identifier ) const;

private:

	char ligand_chain_ = 'X';

	std::string pose_path_;

	std::string main_score_term_;

	core::Size n_runs_ = 0;

	core::Size current_runs_ = 0;

	LigandIdentifier current_ligand_;

	std::string current_ligand_smiles_;

	FragmentLibrary& library_;

	core::pose::PoseOP basic_pose_ = nullptr;

	core::pose::PoseOP best_pose_ = nullptr;

	core::pose::PoseOP working_pose_;

	utility::vector1< moves::MoverOP > mover_;

	std::map< LigandIdentifier, std::map< std::string, core::Real > > score_memory_;

	/// @brief Important! LigendIdentifier here is NOT the same as in the rest of the program. This uses the real names from reactions and reagents
	std::map< std::string, std::map< utility::vector1< std::string >, std::map< std::string, core::Real > > > loaded_score_memory_;

	std::map< LigandIdentifier, std::set< core::Size > > id_memory_;

	bool score_next_ = false;

	core::scoring::ScoreFunctionOP score_function_ = nullptr;

	core::Real base_similarity_penalty_ = 0.2;

	core::Real similarity_penalty_threshold_ = 0.95;

};

// TODO: It might be beneficial to include the pnear term from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0240450 into the calculation
//       It quantifies docking success to a range of [0,1] and could be used as a multiplier, eg: fitness = pnear * lid_root2
utility::vector1< std::string > const score_terms_ { "ligand_interface_delta", "total_REU", "ligand_interface_delta_EFFICIENCY", "lid_root2", "lid_root3" };

}
}

#endif
