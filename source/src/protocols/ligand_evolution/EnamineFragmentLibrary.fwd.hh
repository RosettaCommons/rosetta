// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/EnamineFragmentLibrary.fwd.hh
/// @brief  Forward Class declaration for %EnamineFragmentLibrary
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)


#ifndef INCLUDED_protocols_ligand_evolution_EnamineFragmentLibrary_FWD_HH
#define INCLUDED_protocols_ligand_evolution_EnamineFragmentLibrary_FWD_HH

// Utility headers
#include <core/types.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols{
namespace ligand_evolution{

    class EnamineFragmentLibrary;

    typedef utility::pointer::shared_ptr< EnamineFragmentLibrary > EnamineFragmentLibraryOP;
    typedef utility::pointer::shared_ptr< EnamineFragmentLibrary const > EnamineFragmentLibraryCOP;

    class EnamineReagent;

    typedef utility::pointer::shared_ptr< EnamineReagent > EnamineReagentOP;
    typedef utility::pointer::shared_ptr< EnamineReagent const > EnamineReagentCOP;

    class EnamineReaction;

    typedef utility::pointer::shared_ptr< EnamineReaction > EnamineReactionOP;
    typedef utility::pointer::shared_ptr< EnamineReaction const > EnamineReactionCOP;

    /// @brief LigandIdentifiers are a vector of ints which should encode all information needed by the FragmentLibrary to create a new ligand
    typedef utility::vector1< core::Size > LigandIdentifier;

    /// @brief ReagentSimilarityList is a list of reagent ids combined with a similarity score, mainly used for fragment searches
    typedef utility::vector1< std::pair< core::Size, core::Real > > ReagentSimilarityList;

}
}

#endif
