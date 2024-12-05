// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/Crossover.hh
/// @brief  Class declaration of the %Crossover class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)


#ifndef INCLUDED_protocols_ligand_evolution_Crossover_HH
#define INCLUDED_protocols_ligand_evolution_Crossover_HH


// unit headers
#include <protocols/ligand_evolution/Crossover.fwd.hh>
#include <protocols/ligand_evolution/OffspringFactory.hh>

// package headers
#include <protocols/ligand_evolution/FragmentLibrary.hh>

namespace protocols {
namespace ligand_evolution {

/// @brief Takes two individuals and produces a set amount of offspring using a crossover.
class Crossover : public OffspringFactory {
public:

	explicit Crossover(  FragmentLibrary const& library );

	/// @brief Selects randomly two parents to create offspring until n_offspring is available
	/// @detail One random parent donates the reaction and both donate minimum one reagent. If they come from different reactions, the alien reagent is mapped to the most similar
	///         counterpart.
	utility::vector1< Individual > apply( utility::vector1< Individual > const& parents, core::Size n_offspring ) const override;

	/// @brief Returns the name of this factory
	std::string const& name() const override;

private:

	/// @brief Crosses two individuals into one offspring
	Individual cross( Individual const& reaction_parent, Individual const& other_parent ) const;

private:

	/// @brief Since this function regularly produces offspring and consults a %FragmentLibrary for this, it keeps ownership to one
	FragmentLibrary const& library_;

	std::string name_ = "Crossover";

};

}
}

#endif
