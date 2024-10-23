// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/TournamentSelector.hh
/// @brief  Class declaration of the %TournamentSelector class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

#ifndef INCLUDED_protocols_ligand_evolution_TournamentSelector_HH
#define INCLUDED_protocols_ligand_evolution_TournamentSelector_HH

// unit headers
#include <protocols/ligand_evolution/TournamentSelector.fwd.hh>

// package headers
#include <protocols/ligand_evolution/Selector.hh>

namespace protocols {
namespace ligand_evolution {

/// @brief Selects only the best scored individuals
class TournamentSelector : public Selector {
public:

	TournamentSelector( core::Size tournament_size, core::Real probability );

	/// @brief Returns the indices of the n_select tournament winning individuals in population.
	utility::vector1< Individual > apply( Population& population, core::Size size, bool remove ) const override;

	/// @brief Return the name of this selector
	std::string const& name() const override;

	/// @brief Sets the size for each tournament
	void set_tournament_size( core::Size size );

	/// @brief Sets the probability that the tournament winner accepts
	void set_probability( core::Real probability );

private:

	std::string name_ = "TournamentSelector";

	/// @brief How many individuals will be randomly selected to participate in the tournament
	core::Size tournament_size_ = 6;

	/// @brief At which probability will an individual be selected, starting with the tournament winner
	core::Real probability_ = 1;

};

}
}

#endif
