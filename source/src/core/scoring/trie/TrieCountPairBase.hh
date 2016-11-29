// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/trie/trie_vs_trie.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_trie_TrieCountPairBase_hh
#define INCLUDED_core_scoring_trie_TrieCountPairBase_hh

// Unit Headers
#include <core/scoring/trie/TrieCountPairBase.fwd.hh>

// Package Headers
#include <core/scoring/trie/RotamerTrie.fwd.hh>
#include <core/scoring/trie/RotamerTrieBase.fwd.hh>

// Project Headers
#include <core/scoring/etable/etrie/EtableAtom.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_1.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.fwd.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.fwd.hh>
#include <core/scoring/etable/EtableEnergy.fwd.hh>

#include <core/scoring/lkball/lkbtrie/LKBAtom.fwd.hh>
#include <core/scoring/lkball/lkbtrie/LKBTrieEvaluator.fwd.hh>

#include <core/scoring/elec/electrie/ElecAtom.fwd.hh>
#include <core/scoring/elec/electrie/ElecTrieEvaluator.fwd.hh>

#include <core/scoring/hbonds/hbtrie/HBAtom.fwd.hh>
#include <core/scoring/hbonds/hbtrie/HBCPData.fwd.hh>
#include <core/scoring/hbonds/HBondEnergy.fwd.hh>

#include <core/scoring/methods/MMLJEnergyInter.fwd.hh>
#include <core/scoring/mm/mmtrie/MMEnergyTableAtom.fwd.hh>

#include <core/scoring/vdwaals/VDWTrie.fwd.hh>
#include <core/scoring/vdwaals/VDW_Energy.fwd.hh>

#include <core/types.hh>


// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.fwd.hh>

namespace core {
namespace scoring {
namespace trie {

class TrieCountPairBase : public utility::pointer::ReferenceCount
{
public:

	/// When returning the path distance for two atoms that are separated by more chemical
	/// bonds than any portion of the score function cares, return this as the equivalent
	/// of an infinite chemical separation.
	static const Size INFINITE_SEPARATION = 5;

	virtual ~TrieCountPairBase();
	// There are as many as T*K*M*N^2 different type resolution functions for
	// M different scoring functions with
	// N different count pair data types, spread over
	// K different CountPair classes,
	// for T different Trie-related functions (trie_vs_trie, trie_vs_path, path_vs_path)
	// Not all M*N^2 combinations are reasonable; many of the K*M*N^2 will cause
	// a utility::exit call.
	// Obnoxious, yes, but easier than rewriting the trie_vs_trie algorithm, and
	// worth* the speedup over polymorphism in the innermost loops
	//
	// *speculation.  I will reimplement using polymorphism and compare the speed.
	//
	// Currently, T = 2, M = 1, N = 3, and K = 3

	///---------- TYPE RESOLUTION FUNCTIONS ----------///
	/// What we have here is triple dispatch.  We have to determine the type of
	/// the first rotamer trie, the type of the second rotamer trie, and finally,
	/// the type of the count pair function.  With the type of all three objects known,
	/// then the templated version of the trie-vs-trie method may be invoked (or,
	/// as the compiler sees it, created).  There is no good way to do triple dispatch
	/// (or double dispatch for that matter) in C++.  Polymorphism alows single dispatch,
	/// but nothing further.
	/// The decision I've made was to have the count pair object the last object to resolve
	/// its type.  The count pair object is simpler than the trie objects, and so this ugly
	/// type resolution code polutes what is otherwise a trivial class.  If it were in
	/// the RotamerTrie class, it would be uglier.  As it is, the second trie to identify its
	/// type has to have methods for each of the possible types of the first trie that identifies
	/// its type... its ugly, yes, but not as bad as this code.
	/// Before this design becomes the standard, I have to answer the following two questions:
	/// - Does using templates instead of polymorphism save so much time that its worth this uglyness?
	/// - Does using templates instead of polymorphism increase the size of the executable by an
	///   acceptably small margin?


	//////////////////////////// EtableEnergy /////////////////////////////////////
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrieBase const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);


	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrieBase const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrieBase const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrieBase const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrieBase const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrieBase const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrieBase const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrieBase const & trie2,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	//////////////////////////////////////// EtableEnergy -- analytic evaluation /////////////////////////////////
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrieBase const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrieBase const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrieBase const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrieBase const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrieBase const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrieBase const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrieBase const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrieBase const & trie2,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	//////////////////////////// HBONDS /////////////////////////////////
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & trie1,
		RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & trie2,
		hbonds::HBondEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & trie1,
		RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & trie2,
		hbonds::HBondEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;


	// called when hbonds and etable tries get confused
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & trie1,
		RotamerTrieBase const & trie2,
		hbonds::HBondEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & trie1,
		RotamerTrieBase const & trie2,
		hbonds::HBondEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	/// Hack Elec Energy

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	// called when hack elec tries get confused with other tries
	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrieBase const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrieBase const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrieBase const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrieBase const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrieBase const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrieBase const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrieBase const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< elec::electrie::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrieBase const & trie2,
		elec::electrie::ElecTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	//////////////////////////// LKball //////////////////////////////////////////////
	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	// called when hack elec tries get confused with other tries
	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrieBase const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrieBase const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrieBase const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrieBase const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrieBase const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrieBase const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrieBase const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< lkball::lkbtrie::LKBAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrieBase const & trie2,
		lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	//////////////////////////// MMLJEnergyInter /////////////////////////////////////
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrieBase const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);


	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrieBase const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrieBase const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrieBase const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrieBase const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrieBase const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrieBase const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrieBase const & trie2,
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
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrieBase const & trie2,
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
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrieBase const & trie2,
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
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrieBase const & trie2,
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
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table) = 0;

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrieBase const & trie2,
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
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		RotamerTrieBase const & trie2,
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
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		RotamerTrieBase const & trie2,
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
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		RotamerTrieBase const & trie2,
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
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector) = 0;

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		RotamerTrieBase const & trie2,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


};


} // namespace trie
} // namespace scoring
} // namespace core

#endif
