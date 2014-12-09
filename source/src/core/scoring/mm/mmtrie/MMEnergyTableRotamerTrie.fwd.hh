// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/mmtrie/MMEnergyTableTrie.fwd.hh
/// @brief  Trie class for MMLJEnergy
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_core_scoring_mm_mmtrie_MMEnergyTableRotamerTrie_fwd_hh
#define INCLUDED_core_scoring_mm_mmtrie_MMEnergyTableRotamerTrie_fwd_hh

#include <core/scoring/trie/RotamerTrieBase.fwd.hh>

namespace core {
namespace scoring {
namespace mm {
namespace mmtrie {

typedef trie::RotamerTrieBaseOP MMEnergyTableRotamerTrieOP;
typedef trie::RotamerTrieBaseCOP MMEnergyTableRotamerTrieCOP;

} // namespace mmtrie
} // namespace mm
} // namespace scoring
} // namespace core

#endif
