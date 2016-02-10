// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/trie/trie_vs_trie.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/etable/etrie/TrieCountPairNone.hh>

#include <utility/exit.hh>
#include <iostream>


namespace core {
namespace scoring {
namespace etable {
namespace etrie {

using namespace trie;

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

///////////////////////// EtableEnergy -- analytic evaluation //////////////////////////////////

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table

)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table

)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table

)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table

)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,

	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table

)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,

	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table

)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,

	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table

)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,

	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table

)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,

	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table

)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,

	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table

)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table

)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,

	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table

)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/

)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/

)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const &, //trie2,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


// HBONDS
void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	hbonds::HBondEnergy const & ,
	ObjexxFCL::FArray2D< core::PackerEnergy > & ,
	ObjexxFCL::FArray2D< core::PackerEnergy > & )
{
	utility_exit_with_message( "etable::etrie::TrieCountPairNone::resolve_trie_vs_trie reached with HBondEnergy" );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	hbonds::HBondEnergy const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "etable::etrie::TrieCountPairNone::resolve_trie_vs_path reached with HBondEnergy" );
}


/// Hack Elec E
void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

//////////////////////////// lkball ////////////////////////////////
void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & /*trie2*/,
	lkball::lkbtrie::LKBTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


///////////////////////// MMLJEnergyInter //////////////////////////

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

///////////////////////// VDW Energy //////////////////////////

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_trie(
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	//trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::resolve_trie_vs_path(
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const &, //trie2,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	//trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairNone::print()
{
	std::cout << "CountPairNone" << std::endl;
}


} // namespace etrie
} // namespace etable
} // namespace scoring
} // namespace core
