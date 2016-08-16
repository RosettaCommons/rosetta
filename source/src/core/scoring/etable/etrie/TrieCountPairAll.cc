// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/etrie/TrieCountPairAll.cc
/// @brief  The functions that belong to TrieCountPairAll other than the templated trie-vs-trie
///         and trie-vs-path templated function calls.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/etable/etrie/TrieCountPairAll.hh>

// STL Headers
#include <iostream>


namespace core {
namespace scoring {
namespace etable {
namespace etrie {

using namespace trie;

TrieCountPairAll::~TrieCountPairAll() {}

void
TrieCountPairAll::print()
{
	std::cout << "TrieCountPairAll" << std::endl;
}

} // namespace etrie
} // namespace etable
} // namespace scoring
} // namespace core

