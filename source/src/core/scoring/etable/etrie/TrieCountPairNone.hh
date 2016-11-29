// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/etrie/TrieCountPairNone.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_etable_etrie_TrieCountPairNone_hh
#define INCLUDED_core_scoring_etable_etrie_TrieCountPairNone_hh

// Unit Headers
#include <core/scoring/etable/etrie/TrieCountPairNone.fwd.hh>

// Package Headers
#include <core/scoring/etable/etrie/CountPairData_1_1.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.fwd.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.fwd.hh>

#include <core/scoring/elec/electrie/ElecAtom.fwd.hh>
#include <core/scoring/elec/electrie/ElecTrieEvaluator.fwd.hh>
#include <core/scoring/hbonds/HBondEnergy.fwd.hh>
#include <core/scoring/hbonds/hbtrie/HBCPData.fwd.hh>

#include <core/scoring/lkball/lkbtrie/LKBAtom.fwd.hh>
#include <core/scoring/lkball/lkbtrie/LKBTrieEvaluator.fwd.hh>

#include <core/scoring/trie/TrieCountPairBase.hh>
#include <core/scoring/trie/RotamerTrie.fwd.hh>

#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <core/scoring/etable/etrie/EtableAtom.fwd.hh>

//XRW_B_T1
//#include <core/scoring/etable/CoarseEtableEnergy.fwd.hh>
//XRW_E_T1
#include <core/scoring/methods/MMLJEnergyInter.fwd.hh>

// Project Headers
#include <core/types.hh>

// Objexx FCL Headers
#include <ObjexxFCL/FArray2D.fwd.hh>

namespace core {
namespace scoring {
namespace etable {
namespace etrie {

class TrieCountPairNone : public trie::TrieCountPairBase
{
public:

	///------------ Useful Functions ------------------///
	template < class CPDATA1, class CPDATA2 >
	bool operator () ( CPDATA1 const & /*at1dat*/, CPDATA2 const & /*at2dat*/, Real & /*weight*/, Size & /*path_dist*/ )
	{
		return false;
		//Size const path_dist( at1dat.conn_dist() + at2dat.conn_dist() + 1 );
		//std::cout << "trie_cp1bc3 " << at1dat.conn_dist() << " " << at2dat.conn_dist() << std::endl;
		//return count_pair::CountPairCrossover3::count_at_path_distance( path_dist, weight );
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
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	/////////////////////////////// lkball ///////////////////////////////////
	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
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
