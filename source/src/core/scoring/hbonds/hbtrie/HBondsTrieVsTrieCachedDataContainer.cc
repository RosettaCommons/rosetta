// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/hbonds/hbtrie/HBondsTrieVsTrieCachedDataContainer.cc
/// @brief A class for passing data to the trie-vs-trie calculation for hydrogen bonds, without
/// having to cache it in mutable data in the HBondEnergy method or whatnot.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <core/scoring/hbonds/hbtrie/HBondsTrieVsTrieCachedDataContainer.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "core.scoring.hbonds.hbtrie.HBondsTrieVsTrieCachedDataContainer" );


namespace core {
namespace scoring {
namespace hbonds {
namespace hbtrie {

/// @brief Copy constructor.  Keep default unless deep copying is needed (and in that case,
/// consider using DeepCopyOPs.)
HBondsTrieVsTrieCachedDataContainer::HBondsTrieVsTrieCachedDataContainer( HBondsTrieVsTrieCachedDataContainer const & )=default;

/// @brief Destructor.
HBondsTrieVsTrieCachedDataContainer::~HBondsTrieVsTrieCachedDataContainer(){}

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
core::scoring::trie::TrieVsTrieCachedDataContainerBaseOP
HBondsTrieVsTrieCachedDataContainer::clone() const {
	return utility::pointer::make_shared< HBondsTrieVsTrieCachedDataContainer >( *this );
}

} //hbtrie
} //hbonds
} //scoring
} //core
