// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/selectors/Selector.hh
/// @brief  Class declaration for %Selector
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

#ifndef INCLUDED_protocols_ligand_evolution_Selector_HH
#define INCLUDED_protocols_ligand_evolution_Selector_HH

// unit headers
#include <protocols/ligand_evolution/selectors/Selector.fwd.hh>

// project headers
#include <protocols/ligand_evolution/Population.hh>

// utility headers
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace ligand_evolution {

/// @brief: An abstract class to give an interface for subset selection of %Population
/// @details comming soon.
class Selector : public utility::VirtualBase {
public:

	Selector();
	~Selector();

	/// @brief Applies the specified rules to create a list of indices of selected individuals
	virtual utility::vector1< Individual > apply( Population& population, core::Size size, bool remove ) const = 0;

	/// @brief Return the name of this selector
	virtual std::string const& name() const = 0;

};


}
}

#endif
