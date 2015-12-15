// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/MMLJEnergyInter.hh
/// @brief  molecular mechanics lj energy
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_core_scoring_methods_MMLJEnergyInter_hh
#define INCLUDED_core_scoring_methods_MMLJEnergyInter_hh

// Unit headers
#include <core/scoring/methods/MMLJEnergyInter.fwd.hh>
#include <core/scoring/mm/MMLJEnergyTable.hh>

#include <core/scoring/mm/mmtrie/MMEnergyTableAtom.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>

#include <core/scoring/trie/RotamerTrieBase.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// C++ headers
#include <iostream>

#include <core/scoring/trie/TrieCountPairBase.fwd.hh>
#include <utility/vector1.hh>

#ifdef SERIALIZATION
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace methods {

class MMLJEnergyInter : public ContextIndependentTwoBodyEnergy {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;
public:

	/// ctor
	MMLJEnergyInter();

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	void
	setup_for_minimizing(
		pose::Pose & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & min_map
	) const;


	virtual
	void
	setup_for_packing(
		pose::Pose & pose,
		utility::vector1< bool > const &,
		utility::vector1< bool > const &
	) const;

	// Creates a rotamer trie for the input set of rotamers and stores the trie
	// in the rotamer set.
	virtual
	void
	prepare_rotamers_for_packing(
		pose::Pose const & pose,
		conformation::RotamerSetBase & set ) const;

	// Updates the cached rotamer trie for a residue if it has changed during the course of
	// a repacking
	virtual
	void
	update_residue_for_packing( pose::Pose & pose, Size resid ) const;


	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	virtual
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const & /* domain_map*/,
		ScoreFunction const & /*sfxn*/,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & ) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	virtual
	void
	evaluate_rotamer_pair_energies(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
	) const;

	//@brief overrides default rotamer/background energy calculation and uses
	// the trie-vs-trie algorithm instead
	virtual
	void
	evaluate_rotamer_background_energies(
		conformation::RotamerSetBase const & set,
		conformation::Residue const & residue,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		utility::vector1< core::PackerEnergy > & energy_vector
	) const;

	/// @brief MMLJEnergyInter does not have an atomic interation threshold
	virtual
	Distance
	atomic_interaction_cutoff() const;

	/// @brief MMLJEnergyInter is context independent; indicates that no context graphs are required
	virtual
	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const;

	/// @brief required for neighbor list and to be more lke the ETable
	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		Size res1,
		Size res2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn
	) const;

	/// @brief required for neighbor list and to be more lke the ETable
	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn
	) const;

	/// @brief required for neighbor list and to be more lke the ETable
	etable::count_pair::CountPairFunctionOP
	get_intrares_countpair(
		conformation::Residue const & res,
		pose::Pose const & pose,
		ScoreFunction const & sfxn
	) const;

	trie::RotamerTrieBaseOP
	create_rotamer_trie(
		conformation::RotamerSetBase const & rotset,
		pose::Pose const & pose
	) const;

	trie::RotamerTrieBaseOP
	create_rotamer_trie(
		conformation::Residue const & res,
		pose::Pose const & pose
	) const;

	// Hooks for trie algorithms

	///  How close do two heavy atoms have to be such that their hydrogen atoms might interact?
	///  max heavy-to-hydrogen distance ( MAGIC NUMBER!!!! FIX IT ) + atom-pair interaction distance.
	Real
	hydrogen_interaction_cutoff2() const
	{
		return ( hydrogen_interaction_cutoff() )*( hydrogen_interaction_cutoff() );
	}

	Real
	hydrogen_interaction_cutoff() const
	{
		return ( potential_.max_dist() + 2 * 1.468 ); // HR3 has the largest hydrogen radii @ 1.468
	}

	trie::TrieCountPairBaseOP
	get_count_pair_function_trie(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		trie::RotamerTrieBaseCOP trie1,
		trie::RotamerTrieBaseCOP trie2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn
	) const;

	trie::TrieCountPairBaseOP
	get_count_pair_function_trie(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn
	) const;

	inline
	Energy heavyatom_heavyatom_energy(
		core::scoring::mm::mmtrie::MMEnergyTableAtom const & at1,
		core::scoring::mm::mmtrie::MMEnergyTableAtom const & at2,
		DistanceSquared & d2,
		Size & path_dist
	) const
	{
		// This is not computed higher up in the chain!
		d2 = at1.xyz().distance_squared( at2.xyz() );
		Real rep(0), atr(0);
		potential_.score( at1.mm_atom_type(), at2.mm_atom_type(), path_dist, d2, rep, atr );
		return rep + atr;

	}

	inline
	Energy heavyatom_hydrogenatom_energy(
		core::scoring::mm::mmtrie::MMEnergyTableAtom const & at1,
		core::scoring::mm::mmtrie::MMEnergyTableAtom const & at2,
		Size & path_dist
	) const
	{
		Real dist_squared( at1.xyz().distance_squared( at2.xyz() ) );
		Real rep(0), atr(0);
		potential_.score( at1.mm_atom_type(), at2.mm_atom_type(), path_dist, dist_squared, rep, atr );
		return rep + atr;
	}

	inline
	Energy hydrogenatom_heavyatom_energy(
		core::scoring::mm::mmtrie::MMEnergyTableAtom const & at1,
		core::scoring::mm::mmtrie::MMEnergyTableAtom const & at2,
		Size & path_dist
	) const
	{
		Real dist_squared( at1.xyz().distance_squared( at2.xyz() ) );
		Real rep(0), atr(0);
		potential_.score( at1.mm_atom_type(), at2.mm_atom_type(), path_dist, dist_squared, rep, atr );
		return rep + atr;
	}

	inline
	Energy hydrogenatom_hydrogenatom_energy(
		core::scoring::mm::mmtrie::MMEnergyTableAtom const & at1,
		core::scoring::mm::mmtrie::MMEnergyTableAtom const & at2,
		Size & path_dist
	) const
	{
		Real dist_squared( at1.xyz().distance_squared( at2.xyz() ) );
		Real rep(0), atr(0);
		potential_.score( at1.mm_atom_type(), at2.mm_atom_type(), path_dist, dist_squared, rep, atr );
		return rep + atr;
	}


	// stuff for bump check
	virtual
	void
	bump_energy_full(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	virtual
	void
	bump_energy_backbone(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

private:
	core::scoring::mm::MMLJEnergyTable const & potential_;
	virtual
	core::Size version() const;

};

} // namespace methods
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_methods_MMLJEnergyInter )
#endif // SERIALIZATION

#endif // INCLUDED_core_scoring_methods_MMLJEnergyInter_HH
