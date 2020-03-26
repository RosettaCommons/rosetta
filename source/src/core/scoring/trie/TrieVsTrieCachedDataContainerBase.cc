// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/trie/TrieVsTrieCachedDataContainerBase.cc
/// @brief A generic base class for passing data to the trie-vs-trie calculation, without having to cache
/// it in mutable data in an EnergyMethod or whatnot.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <core/scoring/trie/TrieVsTrieCachedDataContainerBase.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "core.scoring.trie.TrieVsTrieCachedDataContainerBase" );


namespace core {
namespace scoring {
namespace trie {

/// @brief Default constructor.
TrieVsTrieCachedDataContainerBase::TrieVsTrieCachedDataContainerBase() = default;

/// @brief Copy constructor.  Keep default unless deep copying is needed (and in that case,
/// consider using DeepCopyOPs.)
TrieVsTrieCachedDataContainerBase::TrieVsTrieCachedDataContainerBase( TrieVsTrieCachedDataContainerBase const & )=default;

/// @brief Destructor.
TrieVsTrieCachedDataContainerBase::~TrieVsTrieCachedDataContainerBase(){}


} //trie
} //scoring
} //core
