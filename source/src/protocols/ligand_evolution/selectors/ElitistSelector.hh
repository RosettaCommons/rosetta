// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/selector/selectors/ElitistSelector.hh
/// @brief  Class declaration of the %ElitistSelector class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

#ifndef INCLUDED_protocols_ligand_evolution_ElitistSelector_HH
#define INCLUDED_protocols_ligand_evolution_ElitistSelector_HH

// unit headers
#include <protocols/ligand_evolution/selectors/ElitistSelector.fwd.hh>
#include <protocols/ligand_evolution/selectors/Selector.hh>

namespace protocols {
namespace ligand_evolution {

/// @brief Selects only the best scored individuals
class ElitistSelector : public Selector {
public:

	/// @brief Returns the indices of the size fittest individuals in population.
	utility::vector1< Individual > apply( Population& population, core::Size size, bool remove ) const override;

	/// @brief Return the name of this selector
	std::string const& name() const override;

private:

	std::string name_ = "ElitistSelector";

};

}
}

#endif
