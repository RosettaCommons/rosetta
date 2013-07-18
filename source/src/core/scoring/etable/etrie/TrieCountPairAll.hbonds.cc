// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/etrie/TrieCountPairAll.hbonds.cc
/// @brief  File which handles the (failed) dispatch to the TrieCountPairAll class'
///         trie-vs-trie and trie-vs-path functions from the HBondEnergy.  These
///         functions should never be called.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/etable/etrie/TrieCountPairAll.hh>

// Package Headers
#include <core/scoring/trie/trie_vs_trie.hh>
#include <core/scoring/trie/trie_vs_path.hh>
//#include <core/scoring/etable/etrie/CountPairData_1_1.hh>
//#include <core/scoring/etable/etrie/CountPairData_1_2.hh>
//#include <core/scoring/etable/etrie/CountPairData_1_3.hh>
//#include <core/scoring/etable/etrie/CountPairDataGeneric.hh>

// Project Headers
//#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/etable/EtableEnergy.hh>
//#include <core/scoring/hackelec/HackElecEnergy.hh>
//#include <core/scoring/methods/MMLJEnergyInter.hh>
//#include <core/scoring/vdwaals/VDW_Energy.hh>

// STL Headers
#include <iostream>

//#include <core/scoring/etable/etrie/EtableAtom.hh>
#include <core/scoring/trie/RotamerTrie.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace etable {
namespace etrie {

using namespace trie;

// HBONDS
void
TrieCountPairAll::resolve_trie_vs_trie(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	hbonds::HBondEnergy const & ,
	ObjexxFCL::FArray2D< core::PackerEnergy > & ,
	ObjexxFCL::FArray2D< core::PackerEnergy > & )
{
	utility_exit_with_message( "etable::etrie::TrieCountPairAll::resolve_trie_vs_trie reached with HBondEnergy" );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	hbonds::HBondEnergy const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "etable::etrie::TrieCountPairAll::resolve_trie_vs_path reached with HBondEnergy" );
}


} // namespace etrie
} // namespace etable
} // namespace scoring
} // namespace core

