// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/FA_ElecEnergyCD.hh
/// @brief  Electrostatic energy with a distance-dependant dielectric
/// @author Hahnbeom Park

#ifndef INCLUDED_core_scoring_elec_FA_ElecEnergyCD_hh
#define INCLUDED_core_scoring_elec_FA_ElecEnergyCD_hh

// Unit Headers
#include <core/scoring/elec/FA_GrpElecEnergy.fwd.hh>
#include <core/scoring/elec/GroupElec.hh>

// Package headers
#include <core/scoring/elec/ElecAtom.hh>
#include <core/scoring/hbonds/HBondSet.hh>

// Project headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/trie/RotamerTrieBase.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>
#include <core/scoring/trie/TrieCountPairBase.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers

#include <cmath>

namespace core {
namespace scoring {
namespace elec {

class FAElecContextData : public basic::datacache::CacheableData {

public:
 	FAElecContextData();
 	~FAElecContextData();

	void initialize( Size const nres );

	basic::datacache::CacheableDataOP clone() const {
		return FAElecContextDataOP( new FAElecContextData( *this ) );
	}

	Real &n( core::Size i ){ return n_[i]; };
	Real get_n( core::Size i ) const { return n_[i]; };
	Vector &dw_dr( core::Size i ){ return dw_dr_[i]; };
	Vector get_dw_dr( core::Size i ) const { return dw_dr_[i]; };
	utility::vector1< Size > &boundary_neighs( core::Size i ) { return boundary_neighs_[i]; };
	utility::vector1< Size > get_boundary_neighs( core::Size i ) const { return boundary_neighs_[i]; };

private:
	utility::vector1< Real > n_;
	utility::vector1< Vector > dw_dr_;
	utility::vector1< utility::vector1< Size > > boundary_neighs_;
};

///
class FA_GrpElecEnergy : public methods::ContextDependentTwoBodyEnergy  {
public:
	typedef methods::ContextDependentTwoBodyEnergy  parent;
public:

	///
	FA_GrpElecEnergy( methods::EnergyMethodOptions const & options );

	///
	FA_GrpElecEnergy( FA_GrpElecEnergy const & src );

  /// @brief Initilize constants.
	void
	initialize();

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/// stashes nblist if use_nblist is true
	virtual
	void
	setup_for_minimizing(
		pose::Pose & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & min_map
	) const;

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const &scfxn ) const;

	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const &scfxn ) const;

	///
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
		pose::Pose const &,
		conformation::RotamerSetBase & ) const;

	// Updates the cached rotamer trie for a residue if it has changed during the course of
	// a repacking
	virtual
	void
	update_residue_for_packing( pose::Pose &, Size ) const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	/// @brief Returns true if we're using neighborlist-autoupdate
	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & pose ) const;

	virtual
	bool
	defines_score_for_residue_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		bool res_moving_wrt_eachother
	) const;

	virtual
	bool
	use_extended_residue_pair_energy_interface() const;

	virtual
	void
	residue_pair_energy_ext(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResPairMinimizationData const & ,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		kinematics::MinimizerMapBase const &,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData &
	) const;

	/// @brief Evaluate the atom derivative f1/f2 vectors for all atoms on rsd1
	/// in response to the atoms on rsd2, and all the atoms on rsd2 as they
	/// in response to the atoms on rsd1.  This method is used with the
	/// MinimizationGraph and when nblist_autoupdate is not in use.
	virtual
	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const &,
		pose::Pose const & pose, // provides context
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	/// @brief Evaluate the derivative vectors for a particular atom in a given
	/// (asymmetric) pose when nblist_autoupdate is being used.  nblist_autoupdate
	/// cannot be used with symmetric poses, in rtmin, or in minpack.
	virtual
	void
	eval_atom_derivative(
		id::AtomID const &,
		pose::Pose const &,
		kinematics::DomainMap const &,
		ScoreFunction const &,
		EnergyMap const &,
		Vector &,
		Vector &
	) const;

	void
	finalize_total_energy(
		pose::Pose & ,
		ScoreFunction const &,
		EnergyMap & 
	) const;


	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
   ) const;

	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const &,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const;

	virtual
	void
	evaluate_rotamer_pair_energies(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		ScoreFunction const &,
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
		ScoreFunction const &,
		EnergyMap const & weights,
		utility::vector1< core::PackerEnergy > & energy_vector
	) const;


	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	/// @brief Interface function for class NeighborList.
	etable::count_pair::CountPairFunctionCOP
	get_intrares_countpair(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &
	) const;

	/// @brief Interface function for class NeighborList.
	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		Size const,
		Size const,
		pose::Pose const &,
		ScoreFunction const &
	) const;

	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2
	) const;


	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

/// Private methods
private:

	void
	set_nres_mono(
		core::pose::Pose const & pose
	) const;

	void
	precalc_context( pose::Pose & pose,
									 FAElecContextDataOP data
									 ) const;


	Real
	eval_n( Real const cendist,
					Real &dn_drij,
					bool const eval_deriv ) const;

	void
	eval_context_derivatives(
													 conformation::Residue const & rsd1,
													 FAElecContextDataCOP data,
													 EnergyMap const &,
													 utility::vector1< DerivVectorPair > & r1_atom_derivs
													 ) const;

	core::Real
	burial_weight( core::Real const nb ) const;

	core::Real
	burial_deriv( core::Real const nb ) const;

	bool
	monomer_test(
		Size irsd,
		Size jrsd
	) const;

protected:

	inline
	etable::coulomb::Coulomb const &
	coulomb() const {return coulomb_; }

	inline
	GroupElec const &
	groupelec() const {return groupelec_; }

private:

	etable::coulomb::Coulomb coulomb_;
	GroupElec groupelec_;

	bool exclude_protein_protein_;
	bool exclude_monomer_;
	bool exclude_DNA_DNA_;
	Real intrares_scale_;
	bool context_dependent_;

	mutable Size nres_monomer_;

	// temporary for derivative stuffs
	mutable utility::vector1< Real > Eres_;

	virtual
	core::Size version() const;

};

} // namespace elec
} // namespace scoring
} // namespace core

#endif
