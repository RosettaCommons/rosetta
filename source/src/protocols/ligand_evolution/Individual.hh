// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/Individual.hh
/// @brief  Class declaration for %Individual
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)


#ifndef INCLUDED_protocols_ligand_evolution_Individual_HH
#define INCLUDED_protocols_ligand_evolution_Individual_HH

// unit headers
#include <protocols/ligand_evolution/Individual.fwd.hh>

// package headers
#include <protocols/ligand_evolution/EnamineFragmentLibrary.fwd.hh>

// project headers
#include <core/types.hh>

// utility headers
#include <utility/vector1.hh>

// C/C++ headers
#include <map>

namespace protocols {
namespace ligand_evolution {

/// @brief The individual holds all information about a single solution in an evolutionary optimization process
class Individual {

public:

	/// @brief Receives a new pose from the fragment library
	Individual( LigandIdentifier const& identifier, utility::vector1< core::Size > const& parent_ids, std::string const& type_of_birth );
	~Individual();

	/// @brief Returns the score for the given string.
	core::Real score( std::string const& name ) const;

	/// @brief Sets the score according to the given string.
	void score( std::string const& name, core::Real score );

	/// @brief Returns the total score used by the evolutionary optimization
	core::Real score() const;

	/// @brief Sets the total score used by the evolutionary optimization
	void score( core::Real score );

	/// @brief Returns true only if the individual was scored before
	bool is_scored() const;

	/// @brief Returns the LigandIdentifier of this solution
	LigandIdentifier const& identifier() const;

	/// @brief Sets the id for this individual. Returns true if a new id was set
	bool id( core::Size id );

	/// @brief Returns the id of this individual
	core::Size id() const;

	/// @brief Read only access to all raw score terms
	std::map< std::string, core::Real > const& score_terms() const;

	/// @brief Returns this individuals parents ids
	utility::vector1< core::Size > const& parents() const;

	/// @brief Returns this individuals type of birth
	std::string const& type_of_birth() const;

private:

	/// @brief The %LigandIdentifier that describes this Individuals ligand
	LigandIdentifier identifier_;

	/// @brief The total score that is actually used for optimization
	core::Real score_ = 0.0;

	/// @brief Stores different score terms with a string identifier
	std::map< std::string, core::Real > score_terms_;

	/// @brief Saves if a score was already calculated
	bool is_scored_ = false;

	/// @brief A unique number to be identified later
	core::Size id_ = 0;

	/// @brief Stores information about the parents of this individual for later analysis
	utility::vector1< core::Size > parent_ids_;

	/// @brief A short descriptor of how this individual was created
	std::string type_of_birth_;

};

}
}

#endif
