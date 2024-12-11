// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/offspring_factory/IdentityFactory.hh
/// @brief  Class declaration of the %IdentityFactory class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

#ifndef INCLUDED_protocols_ligand_evolution_IdentityFactory_HH
#define INCLUDED_protocols_ligand_evolution_IdentityFactory_HH

// unit headers
#include <protocols/ligand_evolution/offspring_factory/IdentityFactory.fwd.hh>
#include <protocols/ligand_evolution/offspring_factory/OffspringFactory.hh>

namespace protocols {
namespace ligand_evolution {

/// @brief Take a single individual, and returns a desired number of point mutated individuals based on the input one
class IdentityFactory : public OffspringFactory {
public:

	/// @brief Simply copies and returns all individuals. N_offspring won't be used here
	utility::vector1< Individual > apply( utility::vector1< Individual > const& parents, core::Size ) const override;

	/// @brief Returns the name of this mutator
	std::string const& name() const override;

private:

	std::string name_ = "IdentityFactory";

};

}
}

#endif
