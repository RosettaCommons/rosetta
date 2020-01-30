// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/hbonds/hbtrie/HBondsTrieVsTrieCachedDataContainer.fwd.hh
/// @brief A class for passing data to the trie-vs-trie calculation for hydrogen bonds, without
/// having to cache it in mutable data in the HBondEnergy method or whatnot.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_core_scoring_hbonds_hbtrie_HBondsTrieVsTrieCachedDataContainer_fwd_hh
#define INCLUDED_core_scoring_hbonds_hbtrie_HBondsTrieVsTrieCachedDataContainer_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace scoring {
namespace hbonds {
namespace hbtrie {

class HBondsTrieVsTrieCachedDataContainer;

using HBondsTrieVsTrieCachedDataContainerOP = utility::pointer::shared_ptr< HBondsTrieVsTrieCachedDataContainer >;
using HBondsTrieVsTrieCachedDataContainerCOP = utility::pointer::shared_ptr< HBondsTrieVsTrieCachedDataContainer const >;

} //hbtrie
} //hbonds
} //scoring
} //core

#endif //INCLUDED_core_scoring_hbonds_hbtrie_HBondsTrieVsTrieCachedDataContainer_fwd_hh
