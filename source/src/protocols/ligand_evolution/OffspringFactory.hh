// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/OffspringFactory.hh
/// @brief  Class declaration for %OffspringFactory
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

#ifndef INCLUDED_protocols_ligand_evolution_OffspringFactory_HH
#define INCLUDED_protocols_ligand_evolution_OffspringFactory_HH

// unit headers
#include <protocols/ligand_evolution/OffspringFactory.fwd.hh>

// package headers
#include <protocols/ligand_evolution/Individual.hh>

// utility headers
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

namespace protocols{
namespace ligand_evolution{

    /// @brief: An abstract class to give an interface for producing offspring from individuals
    /// @details comming soon.
    class OffspringFactory : public utility::VirtualBase {
    public:

        /// @brief Virtual function to wrap offspring production. Allows for input of any number of individuals and returns n offsprings
        virtual utility::vector1< Individual > apply( utility::vector1< Individual > const& parents, core::Size n_offspring ) const = 0;

        /// @brief Returns the name of this factory
        virtual std::string const& name() const = 0;

    };


}
}

#endif
