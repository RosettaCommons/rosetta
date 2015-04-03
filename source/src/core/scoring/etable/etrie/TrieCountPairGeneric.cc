// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/etrie/TrieCountPairGeneric.hh
/// @brief  Implementation of count pair class to go along with the CountPairGeneric atom data.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/etable/etrie/TrieCountPairGeneric.hh>

// Package Headers
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.hh>
#include <core/scoring/elec/FA_ElecEnergy.hh>
#include <core/scoring/methods/MMLJEnergyInter.hh>
#include <core/scoring/vdwaals/VDW_Energy.hh>

#include <core/scoring/trie/trie_vs_trie.hh>
#include <core/scoring/trie/trie_vs_path.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/conformation/PseudoBond.hh>

// STL Headers
#include <iostream>

#include <core/scoring/etable/etrie/EtableAtom.hh>
#include <core/scoring/trie/RotamerTrie.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace etable {
namespace etrie {

using namespace trie;

/// Grab bond and pseudobond information out of the two residues at construction time.
TrieCountPairGeneric::TrieCountPairGeneric(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	Size res1_cpdat_for_res2,
	Size res2_cpdat_for_res1
) :
	hard_crossover_( false ),
	crossover_( 4 ),
	res1_cpdat_( res1_cpdat_for_res2 ),
	res2_cpdat_( res2_cpdat_for_res1 )
{
	//std::cout << "Constructed Trie Count Pair Generic " << std::endl;
	Size nconnections( 0 );
	Size nbonded( 0 );
	if ( res1.is_bonded( res2.seqpos() )) {
		nbonded = nconnections = res1.connections_to_residue( res2 ).size();
	}
	if ( res1.is_pseudo_bonded( res2.seqpos() )) {
		nconnections += res1.get_pseudobonds_to_residue( res2.seqpos() )->size();
	}

	connection_gaps_.resize( nconnections );
	for ( Size ii = 1; ii <= nbonded; ++ii ) {
		connection_gaps_[ ii ] = 1; // Exactly one bond separates two atoms that are chemically bonded.
	}

	if ( nbonded < nconnections ) {
		using namespace conformation;
		Size count = 0;
		PseudoBondCollectionCOP pbs = res1.get_pseudobonds_to_residue( res2.seqpos() );
		for ( PseudoBondCollection::PBIter pbiter = pbs->iter_begin(), pbiter_end = pbs->iter_end();
				pbiter != pbiter_end; ++pbiter ) {
			++count;
			connection_gaps_[ nbonded + count ] = pbiter->nbonds();
		}
	}
}

void
TrieCountPairGeneric::crossover( Size setting )
{
	crossover_ = setting;
}

void
TrieCountPairGeneric::hard_crossover( bool setting )
{
	hard_crossover_ = setting;
}


///---------- TYPE RESOLUTION FUNCTIONS ----------///

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::TableLookupEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


///////////////////////// EtableEnergy -- analytic evaluation //////////////////////////////////

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/

)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/

)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &  /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(

	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(

	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(

	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & /*trie1*/,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & /*trie2*/,
	etable::AnalyticEtableEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

////////////////////////////////// HBONDS //////////////////////////////////
void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	hbonds::HBondEnergy const & ,
	ObjexxFCL::FArray2D< core::PackerEnergy > & ,
	ObjexxFCL::FArray2D< core::PackerEnergy > & )
{
	utility_exit_with_message( "etable::etrie::TrieCountPairGeneric::resolve_trie_vs_trie reached with HBondEnergy" );
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	hbonds::HBondEnergy const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "etable::etrie::TrieCountPairGeneric::resolve_trie_vs_path reached with HBondEnergy" );
}


/// Hack Elec E
void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}


void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}


void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}


void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}


void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /* temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie2,
	elec::FA_ElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & /*trie2*/,
	elec::FA_ElecEnergy const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie2,
	elec::FA_ElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


///////////////////////// MMLJEnergyInter //////////////////////////

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > &, //pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & //temp_table
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_trie(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const &, //trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const &, //trie2,
	methods::MMLJEnergyInter const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void
TrieCountPairGeneric::resolve_trie_vs_path(
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


/////////////////////////////// VDW Energy ////////////////////////////

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/
)
{
	utility_exit_with_message("Dispatch error. Arrived at TrieCountPairGeneric with incompatible count pair data!");
}

void TrieCountPairGeneric::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & trie2,
	vdwaals::VDWTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairGeneric::print()
{
	std::cout << "TrieCountPairGeneric" << std::endl;
}


} // namespace etrie
} // namespace etable
} // namespace scoring
} // namespace core
