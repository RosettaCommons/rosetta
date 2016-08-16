// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/vdwaals/VDWTrie.fwd.hh
/// @brief  Trie data structure for the low-resolution (centroid) repulsive energy
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_vdwaals_VDWTrie_FWD_HH
#define INCLUDED_core_scoring_vdwaals_VDWTrie_FWD_HH

#include <core/scoring/trie/RotamerTrieBase.fwd.hh>

namespace core {
namespace scoring {
namespace vdwaals {

class VDWAtom;
class VDWTrieCountPair1B;
class VDWTrieCountPairAll;

typedef trie::RotamerTrieBaseOP VDWRotamerTrieOP;
typedef trie::RotamerTrieBaseCOP VDWRotamerTrieCOP;

}
}
}

#endif
