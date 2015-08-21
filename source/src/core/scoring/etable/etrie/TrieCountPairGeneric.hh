// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/etrie/TrieCountPairGeneric.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_etable_etrie_TrieCountPairGeneric_hh
#define INCLUDED_core_scoring_etable_etrie_TrieCountPairGeneric_hh

// Unit Headers
#include <core/scoring/etable/etrie/TrieCountPairGeneric.fwd.hh>

// Package Headers
#include <core/scoring/etable/etrie/CountPairData_1_1.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.fwd.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.hh>

#include <core/scoring/trie/TrieCountPairBase.hh>
#include <core/scoring/trie/RotamerTrie.fwd.hh>

#include <core/scoring/etable/Etable.fwd.hh>
//XRW_B_T1
//#include <core/scoring/etable/CoarseEtableEnergy.fwd.hh>
//XRW_E_T1
#include <core/scoring/elec/FA_ElecEnergy.fwd.hh>
#include <core/scoring/methods/MMLJEnergyInter.fwd.hh>

#include <core/scoring/hbonds/HBondEnergy.fwd.hh>
#include <core/scoring/hbonds/hbtrie/HBCPData.fwd.hh>

// Project Headers
#include <core/types.hh>

// Objexx FCL Headers
#include <ObjexxFCL/FArray2D.fwd.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace etable {
namespace etrie {

class TrieCountPairGeneric : public trie::TrieCountPairBase
{
private:
	/// The number of bonds separating the two residues; for bonds, this is 1, for pseudobonds,
	/// its the number of intervening bonds between the two residue-connection atoms through
	/// the third party residue(s).
	utility::vector1< Size > connection_gaps_;

	// "false" means that atom interactions at the crossover separation should be downweighted but counted.
	// "true" means that atom interactions at the crossover separation are counted at full strength.
	bool hard_crossover_;

	Size crossover_;
	Size const res1_cpdat_;
	Size const res2_cpdat_;


public:

	/// Grab bond and pseudobond information out of the two residues at construction time.
	TrieCountPairGeneric(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		Size res1_cpdat_for_res2,
		Size res2_cpdat_for_res1
	);

	void crossover( Size setting );
	void hard_crossover( bool setting );

	///------------ Useful Functions ------------------///
	bool operator () (
		CountPairDataGeneric const & at1dat,
		CountPairDataGeneric const & at2dat,
		Real & weight,
		Size & path_dist)
	{
		debug_assert( connection_gaps_.size() == at1dat.conn_dat( res1_cpdat_ ).size() );
		debug_assert( connection_gaps_.size() == at2dat.conn_dat( res2_cpdat_ ).size() );
		Size minimum_distance = INFINITE_SEPARATION;
		for ( Size ii = 1; ii <= connection_gaps_.size(); ++ii ) {
			Size iidist = at1dat.conn_dat( res1_cpdat_ )[ ii ] + at2dat.conn_dat( res2_cpdat_ )[ ii ] + connection_gaps_[ ii ];
			if ( iidist < minimum_distance ) minimum_distance = iidist;
		}
		path_dist = minimum_distance;
		if ( !hard_crossover_ && minimum_distance == crossover_ ) weight = count_pair::CountPairFunction::cp_half;
		return minimum_distance >= crossover_;
	}

	static
	void
	print();


	///---------- TYPE RESOLUTION FUNCTIONS ----------///
	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	//////////////////////////////////// EtableEnergy -- analytic evaluation /////////////////////////////////

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	// HBONDS
	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & trie1,
		trie::RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & trie2,
		hbonds::HBondEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);


	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & trie1,
		trie::RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & trie2,
		hbonds::HBondEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	/// Hack Elec Energy


	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	/////////////////////////// MMLJEnergyInter //////////////////////////////
	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	/////////////////////////// VDW_Energy //////////////////////////////
	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

};


} // namespace etrie
} // namespace etable
} // namespace scoring
} // namespace core

#endif
