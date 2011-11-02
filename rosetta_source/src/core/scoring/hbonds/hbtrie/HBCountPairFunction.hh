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
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);


	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);


	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::EtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);



	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::EtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);



//XRW_B_T1
/*
//////////////////////////////////// CoarseEtableEnergy /////////////////////////////////

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);


	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);



	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);


	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);


	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		etable::CoarseEtableEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

*/
//XRW_E_T1

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
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_trie(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_1 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_2 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairData_1_3 > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector);

	virtual
	void
	resolve_trie_vs_path(
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie1,
		trie::RotamerTrie< hackelec::ElecAtom, etable::etrie::CountPairDataGeneric > const & trie2,
		hackelec::HackElecEnergy const & sfxn,
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

};


} // namespace hbtrie
} // namespace hbonds
} // namespace scoring
} // namespace core

#endif
