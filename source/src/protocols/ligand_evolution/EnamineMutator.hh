// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/EnamineMutator.hh
/// @brief  Class declaration of the %EnamineMutator class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

#ifndef INCLUDED_protocols_ligand_evolution_EnamineMutator_HH
#define INCLUDED_protocols_ligand_evolution_EnamineMutator_HH

// unit headers
#include <protocols/ligand_evolution/EnamineMutator.fwd.hh>
#include <protocols/ligand_evolution/OffspringFactory.hh>

// package headers
#include <protocols/ligand_evolution/EnamineFragmentLibrary.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

namespace protocols{
namespace ligand_evolution{

    /// @brief Take a single individual, and returns a desired number of point mutated individuals based on the input one
    class EnamineMutator : public OffspringFactory {
    public:

        EnamineMutator( EnamineFragmentLibrary const& library, utility::vector1< core::Real > const& weights, core::Real min_similarity, core::Real max_similarity );

        /// @brief Accepts a list of parents and creates n_offspring mutants
        /// @detail Based on the internal values this mutation process changes either one of the reagents or the reaction and maps offspring to the closest ligand
        utility::vector1< Individual > apply( utility::vector1< Individual > const& parents, core::Size n_offspring ) const override;

        /// @brief Returns the name of this mutator
        std::string const& name() const override;

        /// @brief Sets the minimum similarity for mutations
        void set_min_similarity( core::Real min_similarity );

        /// @brief Sets the maximum similarity for mutations
        void set_max_similarity( core::Real max_similarity );

    private:

        Individual mutate( Individual const& individual ) const;

    private:

        /// @brief If possible this is the minimum similarity at which a reagent point mutation will be selected
        core::Real min_similarity_ = 0.6;

        /// @brief This is the maximum similarity for reagents to be considered as possible mutants
        core::Real max_similarity_ = 0.99;

        /// @brief Since this function regularly produces offspring and consults a %FragmentLibrary for this, it keeps a ownership to one
        EnamineFragmentLibrary const& library_;

        /// @brief Randomly chooses what feature will be mutated based on the given weights
        numeric::random::WeightedSampler reaction_or_reagent_;

        /// @brief A name to identify the type of this object
        std::string name_ = "EnamineMutator";

    };

}
}

#endif
