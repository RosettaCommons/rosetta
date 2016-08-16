// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/geometric_solvation/OccludedHbondSolEnergy.hh
/// @brief  Solvation model based on penalizing potential for Hbonding to solvent
/// @author John Karanicolas


#ifndef INCLUDED_core_scoring_geometric_solvation_OccludedHbondSolEnergy_hh
#define INCLUDED_core_scoring_geometric_solvation_OccludedHbondSolEnergy_hh

#include <core/types.hh>

// Unit Headers
#include <core/scoring/geometric_solvation/OccludedHbondSolEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/geometric_solvation/DatabaseOccSolEne.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


//#include <core/scoring/EnergyMap.hh>

namespace core {
namespace scoring {
namespace geometric_solvation {

extern Vector dummy_deriv_vector_;

class OccludedHbondSolEnergy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;
public:

	OccludedHbondSolEnergy( methods::EnergyMethodOptions const & options, bool const verbose = false );

	OccludedHbondSolEnergy( OccludedHbondSolEnergy const & src );

	virtual
	methods::EnergyMethodOP
	clone() const;

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & ,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & scorefxn,
		EnergyMap & emap
	) const;

	/// @brief Inform inquiring algorithms that this energy method will opt-in to the
	/// residue-pair decomposable derivative evaluation scheme.
	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	virtual
	bool
	defines_score_for_residue_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		bool res_moving_wrt_eachother
	) const;

	virtual
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
	) const;

	virtual
	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const;

	/*void
	deprecated_eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
	) const;*/

	// Note: This could change - see notes in the .cc re double-counting. If it does,
	// eval_atom_derivative has to change too.
	virtual
	bool
	defines_intrares_energy( EnergyMap const & ) const { return true; };

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > &  ) const {};

private:

	Real
	res_res_occ_sol_one_way(
		conformation::Residue const & polar_rsd,
		conformation::Residue const & occ_rsd ) const;

	void
	eval_residue_pair_derivatives_one_way(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;


	void
	get_atom_atom_occ_solvation(
		Size const don_h_atom,
		Size const don_base_atom,
		conformation::Residue const & don_rsd,
		Size const occ_atom,
		conformation::Residue const & occ_rsd,
		Real & energy,
		bool const update_deriv = false,
		Real const occ_sol_fitted_weight = 0.0,
		//bool const update_deriv_base = false,
		//bool const update_deriv_occ = false,
		Vector & f1_base = dummy_deriv_vector_,
		Vector & f2_base = dummy_deriv_vector_,
		Vector & f1_polar = dummy_deriv_vector_,
		Vector & f2_polar = dummy_deriv_vector_,
		Vector & f1_occ = dummy_deriv_vector_,
		Vector & f2_occ = dummy_deriv_vector_
	) const;

	Real
	get_cos_angle(
		Vector const & base_atom_xyz,
		Vector const & polar_atom_xyz,
		Vector const & occluding_atom_xyz ) const;

	bool
	atom_is_donor_h( conformation::Residue const & rsd, Size const atom ) const;

	bool
	atom_is_acceptor( conformation::Residue const & rsd, Size const atom ) const;

	bool
	atom_is_valid_base( conformation::Residue const & rsd, Size const atom ) const;


private:

	// const-ref to scoring database
	DatabaseOccSolEne const & occ_hbond_sol_database_;

	bool const verbose_;
	virtual
	core::Size version() const;

};

} // geometric_solvation
} // scoring
} // core

#endif

