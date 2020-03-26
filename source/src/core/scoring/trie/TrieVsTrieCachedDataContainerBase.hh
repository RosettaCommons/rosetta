// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/trie/TrieVsTrieCachedDataContainerBase.hh
/// @brief A generic base class for passing data to the trie-vs-trie calculation, without having to cache
/// it in mutable data in an EnergyMethod or whatnot.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_core_scoring_trie_TrieVsTrieCachedDataContainerBase_hh
#define INCLUDED_core_scoring_trie_TrieVsTrieCachedDataContainerBase_hh

#include <core/scoring/trie/TrieVsTrieCachedDataContainerBase.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace core {
namespace scoring {
namespace trie {

/// @brief A generic base class for passing data to the trie-vs-trie calculation, without having to cache
/// it in mutable data in an EnergyMethod or whatnot.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class TrieVsTrieCachedDataContainerBase : public utility::VirtualBase {

public:

	/// @brief Default constructor.
	TrieVsTrieCachedDataContainerBase();

	/// @brief Copy constructor.
	TrieVsTrieCachedDataContainerBase(TrieVsTrieCachedDataContainerBase const & src);

	/// @brief Destructor.
	~TrieVsTrieCachedDataContainerBase() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	/// @details Must be implemented by derived classes.  (This is a pure virtual base class.)
	virtual TrieVsTrieCachedDataContainerBaseOP clone() const = 0;

private:

};

} //trie
} //scoring
} //core

#endif //INCLUDED_core_scoring_trie_TrieVsTrieCachedDataContainerBase_hh
