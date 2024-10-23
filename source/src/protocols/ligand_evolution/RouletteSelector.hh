// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/RouletteSelector.hh
/// @brief  Class declaration of the %RouletteSelector class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

#ifndef INCLUDED_protocols_ligand_evolution_RouletteSelector_HH
#define INCLUDED_protocols_ligand_evolution_RouletteSelector_HH

// unit headers
#include <protocols/ligand_evolution/RouletteSelector.fwd.hh>
#include <protocols/ligand_evolution/Selector.hh>

// numeric headers
#include <numeric/random/WeightedReservoirSampler.hh>

namespace protocols {
namespace ligand_evolution {

/// @brief Selects only the best scored individuals
class RouletteSelector : public Selector {
public:

	/// @brief Returns the indices of the size selected individuals from population with a roulette wheel selection
	utility::vector1< Individual > apply( Population& population, core::Size size, bool remove ) const override;

	/// @brief Return the name of this selector
	std::string const& name() const override;

	/// @brief Switch whether positive fitness should be considered or not. If they are considered, this changes the proportions of all other fitness a bit!
	void consider_positive( bool consider_positive );

private:

	std::string name_ = "RouletteSelector";

	/// @brief Switch whether positive fitness should be considered or not. If they are considered, this changes the proportions of all other fitness a bit!
	bool consider_positive_ = true;

};

}
}

#endif
