// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/magnesium/MgEnergy.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_mg_MgEnergy_HH
#define INCLUDED_core_scoring_mg_MgEnergy_HH

// Unit Headers
#include <core/scoring/magnesium/MgKnowledgeBasedPotential.fwd.hh>

// Package headers
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/func/Func.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace scoring {
namespace magnesium {


class MgEnergy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;

public:


	MgEnergy();

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////
	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn ) const;

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
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const;

	virtual
	void
	setup_for_minimizing_for_residue(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		kinematics::MinimizerMapBase const &,
		basic::datacache::BasicDataCache &,
		ResSingleMinimizationData &
	) const;

	virtual
	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		kinematics::MinimizerMapBase const &,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData & pair_data
	) const;

	virtual
	bool
	requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & ) const;

	virtual
	void
	eval_residue_pair_derivatives(
		conformation::Residue const & ires,
		conformation::Residue const & jres,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	// virtual
	// void
	// eval_atom_derivative(
	//  id::AtomID const & atom_id,
	//  pose::Pose const & pose,
	//  kinematics::DomainMap const & domain_map,
	//  ScoreFunction const & scorefxn,
	//  EnergyMap const & weights,
	//  Vector & F1,
	//  Vector & F2
	// ) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return true; }

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

	virtual
	bool
	use_extended_residue_pair_energy_interface() const { return true; }

	virtual
	void
	residue_pair_energy_ext(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

private:

	void
	residue_pair_energy_one_way(
		conformation::Residue const & rsd1, // The other residue
		conformation::Residue const & rsd2, // The Mg(2+)
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	void
	eval_residue_pair(
		conformation::Residue const & ires,
		conformation::Residue const & jres,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		EnergyMap & emap, // fill score values in here.
		EnergyMap const & weights, // for derivs.
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs) const;

	void
	eval_mg_interaction(
		conformation::Residue const & rsd1 /* other residue */,
		Size const atomno1 /*other atomno */,
		conformation::Residue const & rsd2 /* mg residue */,
		pose::Pose const & pose, // provides context
		EnergyMap & emap,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs /* other residue */,
		utility::vector1< DerivVectorPair > & r2_atom_derivs /* mg residue */
	) const;

	void
	eval_mg_residue_pair_derivatives(
		conformation::Residue const & rsd1 /* other residue */,
		conformation::Residue const & rsd2 /* mg residue */,
		pose::Pose const & pose, // provides context
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs /* other residue */,
		utility::vector1< DerivVectorPair > & r2_atom_derivs /* mg residue */
	) const;

	virtual
	core::Size version() const;

	MgKnowledgeBasedPotentialOP mg_lig_knowledge_based_potential_;
	core::Distance const mg_lig_interaction_cutoff_;
	core::Real const v_angle_width_;
	core::Real const v_angle_width2_;
	core::Real const v_angle_baseline;

	core::Real const mg_ref_score_;
	core::Real const hoh_ref_score_;

	// Following are for solvation. "Wild guesses" from database/chemical/fa_standard/atom_properties.txt --
	// don't use those directly so that we can play with them separately.
	core::Distance const mg_lj_radius_;
	core::Distance const mg_lk_lambda_;
	core::Real     const mg_lk_dgfree_;
	core::Real     const lk_inv_lambda2;
	core::Real     const inv_neg2_tms_pi_sqrt_pi;
	core::Real     const mg_lk_coeff;
	bool const compute_mg_sol_for_hydrogens_;

	core::Distance const mg_sol_interaction_cutoff_;
	core::Distance const mg_sol_fade_zone_;
	func::FuncCOP mg_sol_fade_func_;

};


} //magnesium
} //scoring
} //core


#endif
