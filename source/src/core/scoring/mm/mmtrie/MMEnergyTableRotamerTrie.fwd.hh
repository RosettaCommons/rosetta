// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

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
