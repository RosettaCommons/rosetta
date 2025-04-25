// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/FragmentLibrary.hh
/// @brief  Class declaration for %FragmentLibrary
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)


#ifndef INCLUDED_protocols_ligand_evolution_FragmentLibrary_HH
#define INCLUDED_protocols_ligand_evolution_FragmentLibrary_HH

// unit headers
#include <protocols/ligand_evolution/FragmentLibrary.fwd.hh>

// project headers
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <protocols/ligand_evolution/EvolutionOptions.hh>

// utility headers
#include <utility/vector1.fwd.hh>

// numeric headers
#include <numeric/random/WeightedSampler.hh>

// RDKit headers
#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/Fingerprints/MorganFingerprints.h>
#include <rdkit/GraphMol/RDKitBase.h>

// C/C++ headers
#include <map>
#include <set>

namespace protocols {
namespace ligand_evolution {

/// @brief The %FragmentLibrary implements a combinatorial library for reaction and reagent data. Its main task is to hold chemical information and provide new ligands.
class FragmentLibrary {
public:

	FragmentLibrary();
	~FragmentLibrary();
	// deleted these due to rule of three
	FragmentLibrary& operator=( FragmentLibrary const& other ) = delete;
	FragmentLibrary( FragmentLibrary const& other ) = delete;

public:

	void initialize_from_options( EvolutionOptionsOP options, core::Size external_scoring, core::Size rank );

	/// @brief Sets the internally used pose for ligand creation
	void set_pose( core::pose::PoseCOP pose );

	/// @brief generates a smiles representation for the id and calls create_ligand_pose with smiles
	core::pose::PoseOP create_ligand_pose( LigandIdentifier const& id, bool create_rotamers, char ligand_chain ) const;

	/// @brief generates a new residue represented by a given ligand code and adds it to a detached copy of the internally saved pose object
	core::pose::PoseOP create_ligand_pose( std::string const& smiles, bool create_rotamers, char ligand_chain ) const;

	/// @brief Loads smirks reaction file and all reagents in the same folder
	void load_data( std::string const& reaction_file_path, std::string const& reagent_file_path, core::Size rank );

	/// @brief Loads external smiles for later scoring from file
	void load_smiles( std::string const& path_to_data );

	/// @brief Simple function to generate a completely random ligand
	LigandIdentifier random_ligand() const;

	/// @brief Searches for similar reagents within a given reaction and position
	ReagentSimilarityList get_similar_reagents( core::Size reagent_id, core::Size reaction_id, core::Size position ) const;

	/// @brief Returns the number of reactions
	core::Size reactions_size() const;

	/// @brief Returns the total number of reagents
	core::Size reagents_size() const;

	/// @brief Returns the total number of reagents for the given reaction
	core::Size reagents_size( core::Size reaction_index ) const;

	/// @brief Returns the number of reagents at one specific position in the given reaction
	core::Size reagents_size( core::Size reaction_index, core::Size position ) const;

	/// @brief Returns a random reaction index weighted for size of possible molecules
	core::Size random_reaction() const;

	/// @brief Returns a random reaction that is not included in exclude weighted for size of possible molecules
	core::Size random_reaction( std::set< core::Size > const& exclude ) const;

	/// @brief Returns the id of this reaction
	std::string reaction_id( core::Size reaction_id ) const;

	/// @brief Returns the id of this reagent
	std::string reagent_id( core::Size reagent_id ) const;

	/// @brief Runs the reaction specified by the ligand identifier and returns the resulting rdkit molecule
	std::string run_reaction( LigandIdentifier const& identifier ) const;

	/// @brief Generates a new molecule with rdkit reaction and returns its smiles representation
	std::string identifier_to_smiles( LigandIdentifier const& identifier ) const;

	/// @brief Returns the number of saved external smiles that should be scored
	core::Size n_unscored_smiles() const;

	/// @brief Returns the number of maximum positions for all reactions used
	core::Size max_positions() const;

	/// @brief Returns for a given reaction the number of positions used
	core::Size reaction_positions( core::Size reaction_id ) const;

	/// @brief Returns the used reaction index for a given reaction name. Returns 0 if not found.
	core::Size reaction_name_to_index( std::string const& reaction_name ) const;

	/// @brief Returns the index for a given reagent. Returns 0 if not found.
	core::Size reagent_name_to_index( core::Size reaction_index, core::Size position, std::string const& reagent_name ) const;

	/// @brief Calculates a RDKit Morgan fingerprint for a given ligand. Fingerprints are saved and can be retrieved quickly.
	std::shared_ptr< RDKit::SparseIntVect<unsigned int> > calculate_fingerprint( LigandIdentifier const& id );

	/// @brief Returns the Tanimoto Similarity based on RDKit Morgan fingerprints for two given ligands. The fingerprints are either calculated or retrieved if previously calculated
	core::Real similarity( LigandIdentifier const& id1, LigandIdentifier const& id2 );

private:

	/// @brief LigandIdentifier is interpreted as [reaction, reagent1, reagent2]
	core::conformation::ResidueOP create_ligand( std::string const& smiles, bool create_rotamers ) const;

	/// @brief utility function to load all reactions.
	void load_reactions( std::string const& reaction_file_name );

	/// @brief utility function to load all reagents. Returns indices of newly added entries in reagents_
	void load_reagents( std::string const& reagent_file_path );

	/// @brief Tries to generate rotamers and returns how many were generated
	core::Size generate_rotamers( core::chemical::MutableResidueType& new_ligand ) const;

	/// @brief Takes a ligand and adds it to a detached copy of the internally saved pose object
	core::pose::PoseOP create_pose( core::conformation::Residue& ligand, char ligand_chain ) const;

private:

	/// @brief List of all available reagents
	utility::vector1< ReagentOP > reagents_;

	/// @brief List of all available reactions
	utility::vector1< ReactionOP > reactions_;

	/// @brief Stores the original pose to dock all ligands in
	core::pose::PoseCOP pose_;

	/// @brief Weighted sampler for random unbiased selection
	numeric::random::WeightedSampler weighted_sampler_;

	/// @brief Only used in benchmarking runs. Stores external molecules defined by smiles.
	utility::vector1< std::string > smiles_;

	/// @brief Used during data loading to map reaction names to indices whilst loading reagents
	std::map< std::string, core::Size > reaction_name_to_index_;

	/// @brief stores how many reagents are used at max for a reaction to set the correct length for ligand identifiers
	core::Size max_reagents_ = 0;

	/// @brief stores calculated fingerprints for quick references
	std::map< LigandIdentifier, std::shared_ptr< RDKit::SparseIntVect<unsigned int> > > fingerprints_;

};

/// @brief Internal class to handle reagent information more easily
class Reagent {
public:

	Reagent( std::string const& name, std::string const& reagent_smiles );

	~Reagent() = default;

	// deleted these due to rule of three
	Reagent& operator=( Reagent const& other ) = delete;
	Reagent( Reagent const& other ) = delete;

	std::string name() const { return name_; }

	std::shared_ptr< RDKit::SparseIntVect< unsigned int > > fingerprint() { return fingerprint_; }

	RDKit::RWMOL_SPTR mol() { return mol_; }

private:

	std::string name_;
	RDKit::RWMOL_SPTR mol_;
	std::shared_ptr< RDKit::SparseIntVect< unsigned int > > fingerprint_;

};

/// @brief Internal class to handle reaction information more easily
class Reaction {
public:

	Reaction( std::string const& name, std::string const& reaction_smiles, core::Size n_reagents );

	~Reaction() = default;

	// deleted these due to rule of three
	Reaction& operator=( Reaction const& other ) = delete;
	Reaction( Reaction const& other ) = delete;

	/// @brief returns the index of a random reagent in a FragmentLibraries reagents_ usable for the given position in this reaction
	core::Size random_reagent_index( core::Size position ) const;

	/// @brief Calculates and returns how many molecules possible with this reaction.
	core::Size calculate_possible_molecules();

	/// @brief Returns the number positions with this reaction
	core::Size n_positions() const;

	std::string name() const { return name_; }

	utility::vector1< core::Size > const& reagents( core::Size pos ) const { return reagents_[pos]; }

	void add_reagent( core::Size pos, core::Size reagent_idx ) { reagents_[ pos ].emplace_back( reagent_idx ); }

	core::Size possible_molecules() const { return possible_molecules_; }

	std::shared_ptr< RDKit::ChemicalReaction > reac() { return reac_; }

private:

	std::string name_;
	utility::vector1< utility::vector1< core::Size > > reagents_;
	std::shared_ptr< RDKit::ChemicalReaction > reac_;
	core::Size possible_molecules_ = 1;

};

}
}

#endif
