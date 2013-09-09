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

#ifndef INCLUDED_core_scoring_trie_RotamerTrie_hh
#define INCLUDED_core_scoring_trie_RotamerTrie_hh

// Unit Headers
#include <core/scoring/trie/RotamerTrie.fwd.hh>

// Package Headers
#include <core/scoring/trie/RotamerDescriptor.hh>
#include <core/scoring/trie/RotamerTrieBase.hh>
#include <core/scoring/trie/TrieCountPairBase.hh>
#include <core/scoring/etable/etrie/CountPairData_1_1.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.fwd.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.fwd.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.fwd.hh>
#include <core/scoring/hbonds/hbtrie/HBCPData.hh> // we need full header here because we have inline template function with HBCPData as template specifier

#ifdef WIN32 //VC++ needs full class declaration
 #include <core/scoring/etable/etrie/EtableAtom.hh> // WIN32 INCLUDE
 #include <core/scoring/hbonds/hbtrie/HBAtom.hh> // WIN32 INCLUDE
 #include <core/scoring/mm/mmtrie/MMEnergyTableAtom.hh> // WIN32 INCLUDE
 #include <core/scoring/etable/etrie/CountPairData_1_1.hh> // WIN32 INCLUDE
 #include <core/scoring/etable/etrie/CountPairData_1_2.hh> // WIN32 INCLUDE
 #include <core/scoring/etable/etrie/CountPairData_1_3.hh> // WIN32 INCLUDE
 #include <core/scoring/etable/etrie/CountPairDataGeneric.hh> // WIN32 INCLUDE
 #include <core/scoring/elec/ElecAtom.hh> // WIN32 INCLUDE
 #include <core/scoring/hbonds/hbtrie/HBCPData.hh> // WIN32 INCLUDE
 #include <core/scoring/vdwaals/VDWTrie.hh>
#endif
// Project Headers
#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <core/types.hh>

// Utility Headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/exit.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.fwd.hh>

#include <utility/vector1_bool.hh>


// STL Headers
//#include <cmath>

namespace core {
namespace scoring {
namespace trie {

template < class AT, class CPDATA >
class TrieNode
{
public:
	// CTOR
	TrieNode() :
		subtree_intxn_sphere_radius_sqr_( 0.0 ),
		first_atom_in_branch_( false ),
		is_hydrogen_( false ),
		is_term_( false ),
		sibling_( 0 ),
		rotamers_in_subtree_( 0 )
	{}

	TrieNode( AT atom, CPDATA cp_data ) :
		atom_( atom ),
		cp_data_( cp_data ),
		subtree_intxn_sphere_radius_sqr_( 0.0 ),
		first_atom_in_branch_( false ),
		is_hydrogen_( atom_.is_hydrogen() ),
		is_term_( false ),
		sibling_( 0 ),
		rotamers_in_subtree_( 0 )
	{}

	TrieNode( TrieNode< AT, CPDATA > const & other ) :
		atom_( other.atom_ ),
		cp_data_( other.cp_data_ ),
		subtree_intxn_sphere_radius_sqr_( other.subtree_intxn_sphere_radius_sqr_ ),
		first_atom_in_branch_( other.first_atom_in_branch_ ),
		is_hydrogen_( other.is_hydrogen_ ),
		is_term_( other.is_term_ ),
		sibling_( other.sibling_ ),
		rotamers_in_subtree_( other.rotamers_in_subtree_ )
	{}

	//DSTOR
	~TrieNode() {}

	/// Properties
	bool first_atom_in_branch() const { return first_atom_in_branch_; }
	bool has_sibling() const { return sibling_ != 0; }
	bool is_hydrogen() const { return is_hydrogen_; }
	bool is_rotamer_terminal() const { return is_term_; }

	/// Accessors
	Size num_rotamers_in_subtree() const { return rotamers_in_subtree_; }
	CPDATA const & cp_data() const { return cp_data_; }
	AT const & atom() const { return atom_; }
	Size sibling() const { return sibling_; }
	DistanceSquared subtree_interaction_sphere_square_radius() const {
		return subtree_intxn_sphere_radius_sqr_;
	}

	///Setters
	void first_atom_in_branch( bool setting ) { first_atom_in_branch_ = setting; }
	void sibling( Size setting ) { sibling_ = setting; }
	void is_hydrogen( bool setting ) { is_hydrogen_ = setting; }
	void is_rotamer_terminal( bool setting ) { is_term_ = setting; }
	void subtree_interaction_sphere_square_radius( DistanceSquared setting )
	{
		subtree_intxn_sphere_radius_sqr_ = setting;
	}
	void num_rotamers_in_subtree( Size setting ) { rotamers_in_subtree_ = setting; }

	void print() const
	{
		std::cout << atom_ << "; " << cp_data_ << " : sqrad: " << subtree_intxn_sphere_radius_sqr_;
		std::cout << " faib: " << first_atom_in_branch_;
		std::cout << " H: " << is_hydrogen_;
		std::cout << " term: " << is_term_;
		std::cout << " sib: " << sibling_;
		std::cout << " ris: " << rotamers_in_subtree_ << std::endl;
	}

private:
	AT atom_;
	CPDATA cp_data_;
	DistanceSquared subtree_intxn_sphere_radius_sqr_;
	bool first_atom_in_branch_;
	bool is_hydrogen_;
	bool is_term_;
	Size sibling_;
	Size rotamers_in_subtree_;

};

template < class AT, class CPDATA >
class RotamerTrie : public RotamerTrieBase
{
public:

	RotamerTrie(
		typename utility::vector1< RotamerDescriptor< AT, CPDATA> > & rotamers,
		Distance const atomic_interaction_cutoff_distance
	)
	{
		construct_rotamer_trie( rotamers, atomic_interaction_cutoff_distance );
	}

	virtual ~RotamerTrie() {}

public:   /// Type Resolution Functions
	/// Four trie-vs-trie type resolution functions
	virtual
	void
	trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		other.resolve_trie_vs_trie( *this, cp, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}


	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}

	/// This function is called when the etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrieBase const &,
		TrieCountPairBase &,
		etable::TableLookupEvaluator const &,
		ObjexxFCL::FArray2D< core::PackerEnergy > &,
		ObjexxFCL::FArray2D< core::PackerEnergy > &
	) const
	{
		utility_exit();
	}


	/// Four trie-vs-path type resolution functions
	virtual
	void
	trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		other.resolve_trie_vs_path( *this, cp, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & other,
		TrieCountPairBase & cp,
		etable::TableLookupEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}

	/// This function is called when the etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_path(
		RotamerTrieBase const &,
		TrieCountPairBase & ,
		etable::TableLookupEvaluator const & ,
		utility::vector1< core::PackerEnergy > & ,
		utility::vector1< core::PackerEnergy > &
	) const
	{
		utility_exit_with_message("blah2");
	}

	//////////////////////////////////////// CoarseEtableEnergy//////////////////////////////
	virtual
	void
	trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		other.resolve_trie_vs_trie( *this, cp, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}


	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}


	/// This function is called when the coarse etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrieBase const & ,
		TrieCountPairBase & ,
		etable::AnalyticEtableEvaluator const & ,
		ObjexxFCL::FArray2D< core::PackerEnergy > & ,
		ObjexxFCL::FArray2D< core::PackerEnergy > &
	) const
	{
		utility_exit();
	}

	/// Four trie-vs-path type resolution functions
	virtual
	void
	trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		other.resolve_trie_vs_path( *this, cp, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< etable::etrie::EtableAtom, etable::etrie::CountPairDataGeneric > const & other,
		TrieCountPairBase & cp,
		etable::AnalyticEtableEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}


	/// This function is called when the coarse etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_path(
		RotamerTrieBase const & ,
		TrieCountPairBase & ,
		etable::AnalyticEtableEvaluator const & ,
		utility::vector1< core::PackerEnergy > & ,
		utility::vector1< core::PackerEnergy > &
	) const
	{
		utility_exit();
	}

	// HBond Type Resolution Functions
	virtual
	void
	trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		hbonds::HBondEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		other.resolve_trie_vs_trie( *this, cp, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData > const & other,
		TrieCountPairBase & cp,
		hbonds::HBondEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}

	/// This function is called when hbond energy function get mixed up with non-hbond tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrieBase const & ,
		TrieCountPairBase & ,
		hbonds::HBondEnergy const & ,
		ObjexxFCL::FArray2D< core::PackerEnergy > & ,
		ObjexxFCL::FArray2D< core::PackerEnergy > &
	) const
	{
		utility_exit_with_message( "Type resolution error in trie vs trie.  Unsupported mixing of HBondEnergy function with non-hbond trie.");
	}

	virtual
	void
	trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		hbonds::HBondEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		other.resolve_trie_vs_path( *this, cp, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData > const & other,
		TrieCountPairBase & cp,
		hbonds::HBondEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}

	/// This function is called when hbond energy function get mixed up with non-hbond tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_path(
		RotamerTrieBase const & ,
		TrieCountPairBase & ,
		hbonds::HBondEnergy const & ,
		utility::vector1< core::PackerEnergy > & ,
		utility::vector1< core::PackerEnergy > &
	) const
	{
		utility_exit_with_message("Type resolution error in trie-vs-path.  Unsupported mixing of HBondEnergy function with non-hbond trie.");
	}

	//////////////// Hack Elec E //////////////////////
	/// Four trie-vs-trie type resolution functions
	virtual
	void
	trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		other.resolve_trie_vs_trie( *this, cp, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}


	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}


	/// This function is called when the etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrieBase const &,
		TrieCountPairBase &,
		elec::FA_ElecEnergy const &,
		ObjexxFCL::FArray2D< core::PackerEnergy > &,
		ObjexxFCL::FArray2D< core::PackerEnergy > &
	) const
	{
		utility_exit();
	}


	/// Four trie-vs-path type resolution functions
	virtual
	void
	trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		other.resolve_trie_vs_path( *this, cp, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< elec::ElecAtom, etable::etrie::CountPairDataGeneric > const & other,
		TrieCountPairBase & cp,
		elec::FA_ElecEnergy const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}


	/// This function is called when the etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_path(
		RotamerTrieBase const &,
		TrieCountPairBase & ,
		elec::FA_ElecEnergy const & ,
		utility::vector1< core::PackerEnergy > & ,
		utility::vector1< core::PackerEnergy > &
	) const
	{
		utility_exit_with_message("blah2");
	}

	/// mm lj inter type resolution functions
	virtual
	void
	trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		other.resolve_trie_vs_trie( *this, cp, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}


	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}

	/// This function is called when the etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrieBase const &,
		TrieCountPairBase &,
		methods::MMLJEnergyInter const &,
		ObjexxFCL::FArray2D< core::PackerEnergy > &,
		ObjexxFCL::FArray2D< core::PackerEnergy > &
	) const
	{
		utility_exit();
	}


	/// Four trie-vs-path type resolution functions
	virtual
	void
	trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		other.resolve_trie_vs_path( *this, cp, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< mm::mmtrie::MMEnergyTableAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		methods::MMLJEnergyInter const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}

	/// This function is called when the etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_path(
		RotamerTrieBase const &,
		TrieCountPairBase & ,
		methods::MMLJEnergyInter const & ,
		utility::vector1< core::PackerEnergy > & ,
		utility::vector1< core::PackerEnergy > &
	) const
	{
		utility_exit_with_message("blah2");
	}

	/// VDW_Energy type resolution functions
	virtual
	void
	trie_vs_trie(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		other.resolve_trie_vs_trie( *this, cp, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}

	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}


	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
		ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
	) const
	{
		cp.resolve_trie_vs_trie( other, *this, sfxn, pair_energy_table, temp_table );
	}

	/// This function is called when the etable energy function get mixed up with non-vdwatom tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_trie(
		RotamerTrieBase const &,
		TrieCountPairBase &,
		vdwaals::VDWTrieEvaluator const &,
		ObjexxFCL::FArray2D< core::PackerEnergy > &,
		ObjexxFCL::FArray2D< core::PackerEnergy > &
	) const
	{
		utility_exit_with_message("Type resolution failure in the trie-vs-trie algorithm for the VDW_Energy function");
	}


	/// Four trie-vs-path type resolution functions
	virtual
	void
	trie_vs_path(
		RotamerTrieBase const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		other.resolve_trie_vs_path( *this, cp, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_1 > const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}

	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_2 > const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}


	virtual
	void
	resolve_trie_vs_path(
		RotamerTrie< vdwaals::VDWAtom, etable::etrie::CountPairData_1_3 > const & other,
		TrieCountPairBase & cp,
		vdwaals::VDWTrieEvaluator const & sfxn,
		utility::vector1< core::PackerEnergy > & pair_energy_vector,
		utility::vector1< core::PackerEnergy > & temp_vector
	) const
	{
		cp.resolve_trie_vs_path( other, *this, sfxn, pair_energy_vector, temp_vector );
	}

	/// This function is called when the etable energy function get mixed up with non-etable tries.
	/// It produces a utility_exit call.
	virtual
	void
	resolve_trie_vs_path(
		RotamerTrieBase const &,
		TrieCountPairBase & ,
		vdwaals::VDWTrieEvaluator const & ,
		utility::vector1< core::PackerEnergy > & ,
		utility::vector1< core::PackerEnergy > &
	) const
	{
		utility_exit_with_message("Type resolution failure in the trie-vs-path algorithm for the VDW_Energy function");
	}


	/// END Type Resolution Functions


public:
	/// Useful Functions
	virtual
	void print() const
	{
		std::cout << "RotamerTrie with parameters:" << std::endl;
		for ( Size ii = 1; ii <= num_total_atoms_; ++ii )
		{
			trie_[ ii ].print();
		}
	}

	/// Accessors
	Size num_heavy_atoms() const { return num_heavyatoms_; }
	Size num_unique_rotamers() const { return num_unique_rotamers_; }
	Size max_branch_depth() const { return max_branch_depth_; }
	Size max_heavyatom_depth() const { return max_heavyatom_depth_; }
	Size max_atom_depth() const { return max_atom_depth_; }

	typename utility::vector1< TrieNode< AT, CPDATA > > const &
	atoms() const
	{
		return trie_;
	}

	utility::vector1< Size > const &
	total_rotamers_2_unique_rotamers() const
	{
		return total_rotamers_2_unique_rotamers_;
	}

private: // Functions

	void
	construct_rotamer_trie(
		typename utility::vector1< RotamerDescriptor< AT, CPDATA > > & rotamers,
		Distance const interaction_distance
	)
	{
		using namespace numeric;

		num_total_rotamers_ = rotamers.size();
		total_rotamers_2_unique_rotamers_.resize( num_total_rotamers_ );
		num_total_atoms_ = 0;
		num_heavyatoms_ = 0;

		max_atoms_per_rotamer_ = 0;
		for ( Size ii = 1; ii <= num_total_rotamers_; ++ii ) {
			if ( max_atoms_per_rotamer_ < rotamers[ ii ].natoms() ) {
				max_atoms_per_rotamer_ = rotamers[ ii ].natoms();
			}
		}
		max_atom_depth_ = max_atoms_per_rotamer_; // max_atom_depth_ is redundant!
		std::sort( rotamers.begin(), rotamers.end() );

		utility::vector1< Size > node_stack( max_atoms_per_rotamer_ );
		// index ii represents the number of atoms rotamer ii has with ii-1
		utility::vector1< Size >  num_shared_atoms( num_total_rotamers_, 1);
		for ( Size jj = 2; jj <= num_total_rotamers_; ++jj )
		{
			num_shared_atoms[jj] = rotamers[jj].count_atoms_in_common(rotamers[jj-1]);
			//std::cout << "num_shared_atoms[" << jj << "]: " << num_shared_atoms[jj] << ", natoms[ " << jj-1 << "]: " << rotamers[ jj -1 ].natoms() << std::endl;
		}

		Size count_num_shared_atoms = 0;
		Size num_atoms_in_trie = rotamers[1].natoms();
		for ( Size jj = 2; jj <= num_total_rotamers_; ++jj )
		{
			num_atoms_in_trie += rotamers[jj].natoms() - num_shared_atoms[jj];
			count_num_shared_atoms += num_shared_atoms[jj];
		}
		trie_.resize( num_atoms_in_trie );
		num_total_atoms_ = num_atoms_in_trie;

		//add the first rotamer
		//cerr << "add the first rotamer" << endl;
		max_heavyatom_depth_ = 0;
		utility::vector1< Size > heavyatoms_at_depth( max_atoms_per_rotamer_ );
		heavyatoms_at_depth[1] = 0;

		Size first_rotamer = 1;
		add_atom_to_trie( 1, rotamers[1].atom(1) );

		trie_[1].first_atom_in_branch( true );
		if (! trie_[1].is_hydrogen()) ++num_heavyatoms_;
		node_stack[1] = 1;

		for ( Size jj = 2; jj <= rotamers[first_rotamer].natoms(); ++jj ) {
			//cerr << "atom : " << ii << endl;
			add_atom_to_trie( jj, rotamers[first_rotamer].atom(jj));
			if (! trie_[jj].is_hydrogen()) ++num_heavyatoms_;
			node_stack[jj] = jj;
		}
		trie_[ rotamers[first_rotamer].natoms() ].is_rotamer_terminal( true );

		Size count_unique_rotamers( 1 ), count_atoms_placed( rotamers[1].natoms() );
		total_rotamers_2_unique_rotamers_[rotamers[1].rotamer_id()] = count_unique_rotamers;
		for ( Size jj = 2; jj <= num_total_rotamers_; ++jj ) {
			//int jj_offset_for_bb = total_rotamer_offset_for_bb_(ii) + jj;
			Size jj_num_shared_with_prev = num_shared_atoms[jj];
			Size jj_num_atoms            = rotamers[jj].natoms();
			Size jj_first_distinguished  = jj_num_shared_with_prev + 1;
			//cerr << "rotamer : " << jj << endl;

			if (jj_num_shared_with_prev == jj_num_atoms)
			{
				//duplicate rotamer

				total_rotamers_2_unique_rotamers_[rotamers[jj].rotamer_id()] = count_unique_rotamers;
				continue;
			}

			++count_atoms_placed;
			add_atom_to_trie( count_atoms_placed, rotamers[ jj ].atom( jj_first_distinguished ) );

			// Only mark as the first atom in a brach if there was a genuine branch point
			if ( jj_num_shared_with_prev < rotamers[ jj-1].natoms() ) {
				trie_[count_atoms_placed].first_atom_in_branch( true );
			}

			if (! trie_[count_atoms_placed].is_hydrogen() ) ++num_heavyatoms_;

			assert(  node_stack[jj_first_distinguished]  <= num_atoms_in_trie
				&& (node_stack[jj_first_distinguished]  > 0 || num_shared_atoms[jj] == rotamers[ jj - 1 ].natoms()) );

			if ( node_stack[ jj_first_distinguished] != 0 ) {
				trie_[ node_stack[jj_first_distinguished] ].sibling(count_atoms_placed);
			}
			node_stack[jj_first_distinguished] = count_atoms_placed;

			for ( Size kk = jj_num_shared_with_prev + 2; kk <= jj_num_atoms; ++kk )
			{	//cerr << "atom : " << kk << endl;
				++count_atoms_placed;
				add_atom_to_trie( count_atoms_placed, rotamers[jj].atom(kk) );

				if (! trie_[count_atoms_placed].is_hydrogen()) ++num_heavyatoms_;
				node_stack[ kk ] = count_atoms_placed;
			}

			++count_unique_rotamers;
			total_rotamers_2_unique_rotamers_[rotamers[jj].rotamer_id() ] = count_unique_rotamers;

			trie_[ count_atoms_placed ].is_rotamer_terminal( true );
		}
		num_unique_rotamers_ = count_unique_rotamers;
		compute_max_branch_depth();
		calculate_num_rotamers_in_subtree();
		calculate_subtree_containing_radii( interaction_distance );

	}

	///wow, this function used to be 100 lines long... count pair was such a beast!
	void
	add_atom_to_trie(
		Size trie_atom_id,
		RotamerDescriptorAtom< AT, CPDATA > const & rdatom
	)
	{
		trie_[ trie_atom_id ] = TrieNode< AT, CPDATA >( rdatom.atom(), rdatom.cp_data() );
		return;
	}



	void compute_max_branch_depth()
	{
		//cerr << "begin RotamerTrie::compute_max_stack_height()" << endl;
		max_heavyatom_depth_ = 0;
		max_branch_depth_ = 2;
		utility::vector1< Size > heavy_depth_stack( max_atoms_per_rotamer_ + 1 ); //safe upper bound on size
		heavy_depth_stack[1] = 0;

		Size stack_top = 2;
		for ( Size ii = 1; ii <= num_total_atoms_; ++ii ) {
			if ( trie_[ii].first_atom_in_branch() ) {
				--stack_top;
				//std::cout << "trie_[ " << ii << "].first_atom_in_branch; stack_top: " << stack_top << std::endl;
			}

			if ( trie_[ii].has_sibling() ) {
				//std::cout << "trie_[ " << ii << "].has_sibling: " << trie_[ii].sibling() << " stack_top: " << stack_top + 1 << std::endl;
				++stack_top;
				heavy_depth_stack[stack_top] = heavy_depth_stack[stack_top-1];
			}

			if ( ! trie_[ii].is_hydrogen() ) ++heavy_depth_stack[stack_top];

			if ( max_branch_depth_ < stack_top ) max_branch_depth_ = stack_top;
			if ( max_heavyatom_depth_ < heavy_depth_stack[stack_top] ) {
				max_heavyatom_depth_ = heavy_depth_stack[stack_top];
			}

		}
		//assert( max_branch_depth_ <= max_atoms_per_rotamer_ );
		if( max_branch_depth_ > max_atoms_per_rotamer_+1 ) {
			utility_exit_with_message("max_branch_depth_ should have triggered index-out-of-bounds assertion!");
		}

		//cerr << "max_branch_depth_: " << max_branch_depth_ << endl;
		return;
	}

	// prerequisit: max_branch_depth_ must have been computed
	void calculate_num_rotamers_in_subtree()
	{

		//cerr << "calculate_subtree_heavyatoms_and_rotamers()" << endl;
		utility::vector1< Size > heavyatom_stack( max_atoms_per_rotamer_, 0 );
		utility::vector1< Size > rotamers_in_subtree_stack( max_atoms_per_rotamer_, 0 );

		utility::vector1< Size > heavy_depth_stack( max_branch_depth_ );
		heavy_depth_stack[2] = 0; // why did I set this to -1 before?
		heavy_depth_stack[1] = 0;
		Size stack_top = 2;

		for ( Size ii = 1; ii <= num_total_atoms_; ++ii )
		{
			if ( trie_[ii].first_atom_in_branch() )
			{	for ( Size jj = heavy_depth_stack[stack_top-1] + 1;
					jj <= heavy_depth_stack[stack_top]; ++jj )
				{
					//apl test
					assert ( heavyatom_stack[jj] <= num_total_atoms_ && heavyatom_stack[jj] > 0 );

					trie_[ heavyatom_stack[ jj ]].num_rotamers_in_subtree(
						rotamers_in_subtree_stack[jj]);
					rotamers_in_subtree_stack[jj] = 0;
				}
				--stack_top;
			}

			if ( trie_[ii].has_sibling() ) {
				++stack_top;
				heavy_depth_stack[ stack_top ] = heavy_depth_stack[ stack_top - 1 ];
			}

			if (! trie_[ii].is_hydrogen() ) {
				++heavy_depth_stack[stack_top];
				heavyatom_stack[ heavy_depth_stack[stack_top] ] = ii;
			}

			if ( trie_[ii].is_rotamer_terminal() )	{
				for ( Size jj = 1; jj <= heavy_depth_stack[stack_top]; ++jj ) {
					++rotamers_in_subtree_stack[jj];
				}
			}
		}
		for ( Size ii = 1; ii <= heavy_depth_stack[ stack_top]; ++ii ) {
			assert ( heavyatom_stack[ii] <= num_total_atoms_ && heavyatom_stack[ii] > 0 );
			trie_[ heavyatom_stack[ii] ].num_rotamers_in_subtree( rotamers_in_subtree_stack[ii] );
		}

		//cerr << "done..." << endl;
		return;
	}

	void
	calculate_subtree_containing_radii(
		Distance const interaction_distance
	)
	{
		using namespace numeric;

		//cerr << "calculate_subtree_containing_radii()" << endl;
		utility::vector1< Size > heavyatom_stack( max_atoms_per_rotamer_, 0 );

		utility::vector1< Vector > heavyatom_centers( max_atoms_per_rotamer_ );

		utility::vector1< Size > heavy_depth_stack( max_branch_depth_ + 1 );
		if ( max_branch_depth_ > 0 ) heavy_depth_stack[2] = 0;
		heavy_depth_stack[1] = 0;
		Size stack_top = 2;


		utility::vector1< core::Real > maxd2_in_subtree_stack( max_atoms_per_rotamer_, 0.0 );

		for ( Size ii = 1; ii <= num_total_atoms_; ++ii )
		{
			if ( trie_[ii].first_atom_in_branch() ) {
				for ( Size jj = heavy_depth_stack[stack_top-1] + 1;
					jj <= heavy_depth_stack[stack_top]; ++jj ) {
					core::Real subtree_plus_interaction_diameter =
						std::sqrt(maxd2_in_subtree_stack[ jj ]) + interaction_distance;

					trie_[ heavyatom_stack[ jj ]].subtree_interaction_sphere_square_radius(
						subtree_plus_interaction_diameter * subtree_plus_interaction_diameter );
					maxd2_in_subtree_stack[jj] = 0;
				}
				--stack_top;
			}

			if ( trie_[ii].has_sibling() ) {
				++stack_top;
				heavy_depth_stack[ stack_top ] = heavy_depth_stack[ stack_top - 1 ];
			}

			if ( ! trie_[ii].is_hydrogen() ) {
				++heavy_depth_stack[stack_top];
				heavyatom_stack[ heavy_depth_stack[stack_top] ] = ii;
				heavyatom_centers[ heavy_depth_stack[stack_top]] = trie_[ii].atom().xyz();

				for ( Size jj = 1; jj <= heavy_depth_stack[stack_top] - 1; ++jj ) {
					core::Real const d2 =
						heavyatom_centers[jj].distance_squared(
						heavyatom_centers[heavy_depth_stack[stack_top]] );
					if (d2 > maxd2_in_subtree_stack[jj]) {
						maxd2_in_subtree_stack[jj] = d2;
					}
				}
			}

		}
		for ( Size ii = 1; ii <= heavy_depth_stack[ stack_top]; ++ii ) {
			core::Real subtree_plus_interaction_diameter =
				std::sqrt(maxd2_in_subtree_stack[ ii ]) + interaction_distance;

			trie_[ heavyatom_stack[ii] ].subtree_interaction_sphere_square_radius(
				subtree_plus_interaction_diameter * subtree_plus_interaction_diameter );
		}

		return;
	}


private: // DATA
	typename utility::vector1< TrieNode< AT, CPDATA > > trie_;
	Size num_total_atoms_;
	Size num_heavyatoms_;
	Size max_atoms_per_rotamer_;

	Size num_unique_rotamers_;
	Size num_total_rotamers_;
	utility::vector1< Size > total_rotamers_2_unique_rotamers_;

	Size max_branch_depth_; //rename this max branch depth
	Size max_heavyatom_depth_;
	Size max_atom_depth_;

};


} // namespace trie
} // namespace scoring
} // namespace core

#endif
