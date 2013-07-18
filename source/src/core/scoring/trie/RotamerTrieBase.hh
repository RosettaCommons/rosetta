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

#ifndef INCLUDED_core_scoring_trie_RotamerTrieBase_hh
#define INCLUDED_core_scoring_trie_RotamerTrieBase_hh

// Unit Headers
#include <core/scoring/trie/RotamerTrieBase.fwd.hh>

#include <core/scoring/trie/RotamerTrie.fwd.hh>
#include <core/scoring/trie/TrieCountPairBase.fwd.hh>

/// Package Headers
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/etable/etrie/EtableAtom.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_1.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.fwd.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.fwd.hh>

#include <core/scoring/elec/ElecAtom.fwd.hh>
#include <core/scoring/elec/FA_ElecEnergy.fwd.hh>

#include <core/scoring/hbonds/hbtrie/HBAtom.fwd.hh>
#include <core/scoring/hbonds/HBondEnergy.fwd.hh>
#include <core/scoring/hbonds/hbtrie/HBCPData.fwd.hh>

#include <core/scoring/mm/mmtrie/MMEnergyTableAtom.fwd.hh>
#include <core/scoring/methods/MMLJEnergyInter.fwd.hh>

#include <core/scoring/vdwaals/VDW_Energy.fwd.hh>
#include <core/scoring/vdwaals/VDWTrie.fwd.hh>

// Project Headers
#include <core/conformation/AbstractRotamerTrie.hh>

#include <core/scoring/etable/EtableEnergy.fwd.hh>
//XRW_B_T1
//#include <core/scoring/etable/CoarseEtableEnergy.fwd.hh>
//XRW_E_T1
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.fwd.hh>

// STL Headers
#include <map>

//Auto Headers
namespace core {
namespace scoring {
namespace trie {

class RotamerTrieBase : public conformation::AbstractRotamerTrie
{
public:
	virtual ~RotamerTrieBase() {}

	/// Useful Functions
	virtual
	void print() const = 0;

	//virtual
	//void
	//point_count_pair_data_at_residue( Size other_residue ) const = 0;

	/// The trie contains mutable data to record which peice of count-pair data it is
	/// to use in the upcoming/current trie-vs-trie application.  This data could be held
	/// externally to the trie, and passed in as a parameter to the tvt call, but I believe
	/// the code will be faster if there is no ambiguity over where to look for count pair
	/// information.
	///
	/// Depricated
	//virtual
	//void
	//set_count_pair_data_to_use(
	//	Size connection_id
	//) const = 0;

public:
	void
	set_resid_2_connection_entry( Size resid, Size connid ) {
		resid_2_connid_map_[ resid ] = connid - 1; // index from 0 -- adjust indices here.
	}

	bool
	connection_exists_to_residue( Size resid ) const {
		return resid_2_connid_map_.find( resid ) != resid_2_connid_map_.end();
	}

	Size
	connection_id_for_residue( Size resid ) const {
		/*
		std::cout << "map contents: " << std::endl;
		for ( std::map< Size, Size >::const_iterator iter = resid_2_connid_map_.begin(),
			iter_end = resid_2_connid_map_.end(); iter != iter_end; ++iter ) {
			std::cout << iter->first << " " << iter->second << std::endl;
		}
		*/
		return resid_2_connid_map_.find( resid )->second;
	}


	Size
	get_count_pair_data_for_residue( Size other_residue ) const
	{
		if ( ! connection_exists_to_residue( other_residue ) ) return 0;

		Size connid = connection_id_for_residue( other_residue );
		//std::cout << "other_residue: " << other_residue << " connid: " << connid << std::endl;
		return connid;
	}

private:

	std::map< Size, Size > resid_2_connid_map_;

public:
	// Type resolution Function
	virtual
	void
	trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;


	/// This function is called when the etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	/// Four trie-vs-path type resolution functions
	virtual
	void
	trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;



	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;

	/// This function is called when the etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	////////////////////////////////////////////////////////////////////////////////////////////
	//////////// the same methods again --- overloaded to accept a AnalyticEnergyEvaluator /////
	////////////////////////////////////////////////////////////////////////////////////////////

	virtual
	void
	trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;


	/// This function is called when the coarse etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	/// Four trie-vs-path type resolution functions
	virtual
	void
	trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;



	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	/// This function is called when the coarse etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;

	//// The same methods, again, for the HBondEnergies

	virtual
	void
	trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		hbonds::HBondEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;


	/// This function is called when hbond energy function get mixed up with non-hbond tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		hbonds::HBondEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData > const & other,
		TrieCountPairBase & cp,
		hbonds::HBondEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		hbonds::HBondEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;

	/// This function is called when hbond energy function gets mixed up with non-hbond tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		hbonds::HBondEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData > const & other,
		TrieCountPairBase & cp,
		hbonds::HBondEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;

/// The same methods again, for Hack Elec E
	virtual
	void
	trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;


	/// This function is called when the etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	/// Four trie-vs-path type resolution functions
	virtual
	void
	trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;



	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;



	/// This function is called when non-elec tries get mixed up with the FA_ElecEnergy function
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;

///////////////////////////////////////////////////////////////////////////////////////
//////////// the same methods again --- overloaded for MMLJEnergyTable ////////////////
///////////////////////////////////////////////////////////////////////////////////////

	virtual
	void
	trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	/// This function is called when the mm lj inter energy function get mixed up with other tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	/// Four trie-vs-path type resolution functions
	virtual
	void
	trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;



	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;

	/// This function is called when the mm lj inter energy function get mixed up with other tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;

///////////////////////////////////////////////////////////////////////////////////////
//////////// the same methods again --- overloaded for VDW_Energy /////////////////////
///////////////////////////////////////////////////////////////////////////////////////

	virtual
	void
	trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< core::scoring::vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	/// This function is called when the mm lj inter energy function get mixed up with other tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const = 0;

	/// Four trie-vs-path type resolution functions
	virtual
	void
	trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;



	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< core::scoring::vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;

	/// This function is called when the mm lj inter energy function get mixed up with other tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const = 0;


};


} // namespace trie
} // namespace scoring
} // namespace core

#endif
