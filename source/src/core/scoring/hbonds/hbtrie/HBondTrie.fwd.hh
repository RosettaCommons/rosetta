// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/EtableTrie.fwd.hh
/// @brief  Trie class for etable typedef
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_hbonds_hbtrie_HBondTrie_fwd_hh
#define INCLUDED_core_scoring_hbonds_hbtrie_HBondTrie_fwd_hh

#include <core/scoring/trie/RotamerTrieBase.fwd.hh>

#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace scoring {
namespace hbonds {
namespace hbtrie {

typedef trie::RotamerTrieBaseOP HBondRotamerTrieOP;
typedef trie::RotamerTrieBaseCOP HBondRotamerTrieCOP;

} // hbtrie
} // hbonds
} // scoring
} // core

#endif
