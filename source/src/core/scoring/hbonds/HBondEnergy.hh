// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/HBondEnergy.hh
/// @brief  Hydrogen bond energy method class declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_hbonds_HBondEnergy_hh
#define INCLUDED_core_scoring_hbonds_HBondEnergy_hh

// Unit Headers
#include <core/scoring/hbonds/HBondEnergy.fwd.hh>
//pba

// Package headers
#include <core/scoring/hbonds/hbtrie/HBAtom.hh>
#include <core/scoring/hbonds/hbtrie/HBondTrie.fwd.hh>
#include <core/scoring/hbonds/constants.hh>
//pba

#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/trie/TrieVsTrieCachedDataContainerBase.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <utility/vector1.hh>

#ifdef PYROSETTA
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#endif

#ifdef    SERIALIZATION
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace hbonds {


class HBondEnergy : public methods::ContextDependentTwoBodyEnergy  {
public:
	typedef methods::ContextDependentTwoBodyEnergy  parent;
public:


	HBondEnergy( HBondOptions const & opts );


	HBondEnergy( HBondEnergy const & src );

	~HBondEnergy() override;

	/// clone
	methods::EnergyMethodOP
	clone() const override;


	void
	setup_for_packing(
		pose::Pose & pose,
		utility::vector1< bool > const &,
		utility::vector1< bool > const & ) const override;

	// Creates a rotamer trie for the input set of rotamers and stores the trie
	// in the rotamer set.
	void
	prepare_rotamers_for_packing(
		pose::Pose const & pose,
		conformation::RotamerSetBase & set ) const override;

	// Updates the cached rotamer trie for a residue if it has changed during the course of
	// a repacking
	void
	update_residue_for_packing( pose::Pose & pose, Size resid ) const override;


	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const override;


	/*virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;*/


	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	/// note that this only evaluates sc-sc and sc-bb energies
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const override;

	/// @brief Returns false if two residues are not moving wrt each other; the two parts
	/// of the HBondEnergy function which are non-pairwise-decomposable are held fixed
	/// during minimization -- the neighbor counts, and the bb/bb hbond availability status.
	/// This means that the hbond-energy function can be efficiently evaluated during minimization.
	bool
	defines_score_for_residue_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		bool res_moving_wrt_eachother
	) const override;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override;

	/// @brief Use the extended residue pair energy interface to distinguish between
	/// score function evaluation during minimization from score function evaluation
	/// during regular scoring.
	bool
	use_extended_residue_pair_energy_interface() const override;


	/// @brief Evaluate the energy between a pair of residues during minimization;
	/// during minimization, the bb/bb hbond status is held fixed, so it is possible
	/// to evaluate the bb/bb, bb/sc and sc/sc hydrogen bonds in this function call.
	void
	residue_pair_energy_ext(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;

	/// @brief Setup the bb/bb hbond presence data for a particular residue -- this data
	/// is taken out of the HbondSet in the Pose.
	void
	setup_for_minimizing_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & minmap,
		basic::datacache::BasicDataCache &,
		ResSingleMinimizationData & res_data_cache
	) const override;

	/// @brief Link the bb/bb hbond information in the ResidueSingleMinimizationData
	/// to the ResiduePairMinimizationData.
	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & minmap,
		ResSingleMinimizationData const & res1_data_cache,
		ResSingleMinimizationData const & res2_data_cache,
		ResPairMinimizationData & data_cache
	) const override;

	/// @brief Construct the set of all hydrogen bonds between two residues before
	bool
	requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & pose ) const override;

	/// @brief Do any setup work necessary before evaluating the derivatives for this residue
	/*virtual
	void
	setup_for_derivatives_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const & minsingle_data1,
	ResSingleMinimizationData const & minsingle_data2,
	pose::Pose const & pose,
	ResPairMinimizationData & data_cache
	) const;*/

	/// @brief Retrieves the cached hbond data from the minpair_data object
	/// and calculates the derivative for an atom on rsd1 wrt rsd2.
	/// This method requires that setup_for_derivatives_for_residue_pair
	/// have been called on this residue pair beforehand, and that
	/// the two residues have not changed since that call.
	/*virtual
	void
	eval_atom_derivative_for_residue_pair(
	Size const atom_index,
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const & minsingle_data1,
	ResSingleMinimizationData const & minsingle_data2,
	ResPairMinimizationData const & minpair_data,
	pose::Pose const & pose, // provides context
	kinematics::DomainMap const & domain_map,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
	) const;*/

	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const override;


	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const override;


	//pba
	//virtual
	void
	hbond_derivs_1way(
		EnergyMap const & weights,
		HBondSet const & hbond_set,
		HBondDatabaseCOP database,
		pose::Pose const & pose,
		conformation::Residue const & don_rsd,
		conformation::Residue const & acc_rsd,
		Size const don_nb,
		Size const acc_nb,
		bool const exclude_bsc, /* exclude if acc=bb and don=sc */
		bool const exclude_scb, /* exclude if acc=sc and don=bb */
		Real const ssdep_weight_factor,
		// output
		utility::vector1< DerivVectorPair > & don_atom_derivs,
		utility::vector1< DerivVectorPair > & acc_atom_derivs,
		bool bond_near_wat = false
	) const;


	/// @brief Evaluates the interaction between the backbone of rsd1 and the
	/// backbone of rsd2 and accumulates the unweighted energy.
	void
	backbone_backbone_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;


	/// @brief Evaluates the interaction between the backbone of rsd1 and the
	/// sidechain of rsd2 and accumulates the unweighted energy.
	void
	backbone_sidechain_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;

	/// @brief Evaluates the interaction between the sidechain of rsd1 and the
	/// sidechain of rsd2 and accumulates the unweighted energy.
	void
	sidechain_sidechain_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;


	void
	evaluate_rotamer_pair_energies(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
	) const override;


	//@brief overrides default rotamer/background energy calculation and uses
	// the trie-vs-trie algorithm instead
	void
	evaluate_rotamer_background_energies(
		conformation::RotamerSetBase const & set,
		conformation::Residue const & residue,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		utility::vector1< core::PackerEnergy > & energy_vector
	) const override;


	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const override;

	/// f1 and f2 are zeroed
	/*virtual
	void
	eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
	) const;*/


	Distance
	atomic_interaction_cutoff() const override;


	bool
	divides_backbone_and_sidechain_energetics() const override
	{
		return true;
	}

	Real
	hydrogen_interaction_cutoff2() const;

	/// @brief HBondEnergy is context sensitive
	void indicate_required_context_graphs(
		utility::vector1< bool > & context_graphs_required ) const override;

	bool
	defines_intrares_energy( EnergyMap const & weights ) const override;

	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;

	inline
	Energy heavyatom_heavyatom_energy(
		hbtrie::HBAtom const & at1,
		hbtrie::HBAtom const & at2,
		DistanceSquared & d2,
		Size & /*path_dist*/
	) const
	{
		d2 = at1.xyz().distance_squared( at2.xyz() );
		return 0.0;
	}

	inline
	Energy heavyatom_hydrogenatom_energy(
		hbtrie::HBAtom const & at1, // atom 1 is the heavy atom, the acceptor unless it's a placeholder atom
		hbtrie::HBAtom const & at2, // atom 2 is the hydrogen atom, the donor,
		bool const flipped,
		core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
	) const
	{
		DistanceSquared d2 = at1.xyz().distance_squared( at2.xyz() );
		if ( d2 > MAX_R2 || d2 < MIN_R2 || at1.non_hbonding_atom() ) return 0.0;

		return drawn_out_heavyatom_hydrogenatom_energy( at1, at2, flipped, cached_data );
	}

	Energy
	drawn_out_heavyatom_hydrogenatom_energy(
		hbtrie::HBAtom const & at1, // atom 1 is the heavy atom, the acceptor
		hbtrie::HBAtom const & at2, // atom 2 is the hydrogen atom, the donor
		bool flipped,
		core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
	) const;

	inline
	Energy hydrogenatom_heavyatom_energy(
		hbtrie::HBAtom const & at1,
		hbtrie::HBAtom const & at2,
		Size & /*path_dist*/,
		core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
	) const
	{
		return heavyatom_hydrogenatom_energy( at2, at1, true, cached_data );
	}

	inline
	Energy hydrogenatom_hydrogenatom_energy(
		hbtrie::HBAtom const &,
		hbtrie::HBAtom const &,
		Size & /*path_dist*/,
		core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
	) const
	{
		return 0.0;
	}

	/// @brief Has pairwise energies, but no single-atom ones.
	bool
	has_atomistic_pairwise_energies() const override { return true; }

	void
	atomistic_pair_energy(
		core::Size , // Which atom in residue 1?
		conformation::Residue const &, // Residue 1
		core::Size , // Which atom in residue 2?
		conformation::Residue const &, // Residue 2
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const override;

private:

	hbtrie::HBondRotamerTrieOP
	create_rotamer_trie(
		conformation::RotamerSetBase const & rotset,
		pose::Pose const & pose
	) const;

	hbtrie::HBondRotamerTrieOP
	create_rotamer_trie(
		conformation::Residue const & res,
		pose::Pose const & pose
	) const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	HBondOptionsCOP options_;
	HBondDatabaseCOP database_;

	// 2007 based membrane potentials or new membrane framework initialization
	//
	// (As far as I can tell, this stuff is threadsafe, since these values are only set in calls
	// to setup_for_scoring or setup_for_packing, and are only read after that.  --VKM 22 Nov 2019.)
	mutable Vector normal_;
	mutable Vector center_;
	mutable Real thickness_;
	mutable Real steepness_;
	mutable Real membrane_core_;

	// Used in the "const" evaluate_rotamer_pair_energies and evaluate_rotamer_background_energies
	// methods to keep track of the sfxn weights and the neighbor counts for the two input residues
	// so that the data may be retrieved from within the trie-vs-trie and trie-vs-path calls.
	// mutable EnergyMap weights_;
	// mutable int rotamer_seq_sep_;
	// mutable int res1_;
	// mutable int res2_;
	// mutable int res1_nb_;
	// mutable int res2_nb_;

	core::Size version() const override;

};

} // hbonds
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_hbonds_HBondEnergy )
#endif // SERIALIZATION

#endif

