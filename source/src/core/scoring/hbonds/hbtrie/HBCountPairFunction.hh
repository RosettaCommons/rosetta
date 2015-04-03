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

#ifndef INCLUDED_core_scoring_hbonds_hbtrie_HBCountPairFunction_hh
#define INCLUDED_core_scoring_hbonds_hbtrie_HBCountPairFunction_hh

// Unit Headers
#include <core/scoring/hbonds/hbtrie/HBCountPairFunction.fwd.hh>

// Package Headers
#include <core/scoring/etable/etrie/CountPairData_1_1.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.fwd.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.fwd.hh>
#include <core/scoring/etable/etrie/EtableAtom.fwd.hh>

#include <core/scoring/trie/TrieCountPairBase.hh>
#include <core/scoring/trie/RotamerTrie.fwd.hh>

#include <core/scoring/hbonds/hbtrie/HBCPData.fwd.hh>

#include <core/scoring/hbonds/HBondEnergy.fwd.hh>

// ObjexxFLC Headers
#include <ObjexxFCL/FArray2D.fwd.hh>

namespace core {
namespace scoring {
namespace hbonds {
namespace hbtrie {

class HBCountPairFunction : public trie::TrieCountPairBase
{
public:
	virtual ~HBCountPairFunction();
	/// ------- USEFUL FUNCTIONS -------------///

	/// @brief This is the main function for enforcing the sc/bb hbond exclusion rule.
	/// It works like this: the templated trie-vs-trie function asks this class (via it's operator() method)
	/// whether two atoms should have their interactions counted.  This is answered by the logic in here
	/// that depends on two pieces of data.  1. whether an atom is a sidechain atom, and 2. whether 
	/// the other atom ought to avoid hbonds to sidechain atoms.  This second boolean is true iff
	/// a. the other atom is a backbone atom, b. the other atom is already participating in a bb/bb hbond, and
	/// c. the sc/bb-hydrogen-bond-exclusion rule is being enforced.
	template < class CPDATA1, class CPDATA2 >
	bool operator () ( CPDATA1 at1dat, CPDATA2 at2dat, Real & /*weight*/, Size & /*path_dist*/ )
	{
		return ! ((at1dat.avoid_sc_hbonds() && at2dat.is_sc()) || (at2dat.avoid_sc_hbonds() && at1dat.is_sc()) );
	}

	static
	void
	print();

	///---------- TYPE RESOLUTION FUNCTIONS ----------///
	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);


	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);


	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);


	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	//////////////////////////////////// CoarseEtableEnergy /////////////////////////////////

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);


	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);


	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);


	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	/// HBONDS
	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< HBAtom, HBCPData > const & trie1,
		trie::RotamerTrie< HBAtom, HBCPData > const & trie2,
		HBondEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< HBAtom, HBCPData > const & trie1,
		trie::RotamerTrie< HBAtom, HBCPData > const & trie2,
		HBondEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	/// Hack Elec Energy
	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	/////////////////////////// MMLJEnergyInter //////////////////////////////
	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
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


} // namespace hbtrie
} // namespace hbonds
} // namespace scoring
} // namespace core

#endif
