// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/energy_methods/MgEnergy.hh
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
#include <core/scoring/func/Func.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace energy_methods {


class MgEnergy : public core::scoring::methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentTwoBodyEnergy  parent;

public:


	MgEnergy();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////
	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & scfxn ) const override;

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;


	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &
	) const override;

	void
	setup_for_minimizing_for_residue(
		conformation::Residue const &,
		pose::Pose const &,
		core::scoring::ScoreFunction const &,
		kinematics::MinimizerMapBase const &,
		basic::datacache::BasicDataCache &,
		core::scoring::ResSingleMinimizationData &
	) const override;

	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		kinematics::MinimizerMapBase const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResPairMinimizationData & pair_data
	) const override;

	bool
	requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & ) const override;

	void
	eval_residue_pair_derivatives(
		conformation::Residue const & ires,
		conformation::Residue const & jres,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const override;

	// virtual
	// void
	// eval_atom_derivative(
	//  id::AtomID const & atom_id,
	//  pose::Pose const & pose,
	//  kinematics::DomainMap const & domain_map,
	//  core::scoring::ScoreFunction const & scorefxn,
	//  core::scoring::EnergyMap const & weights,
	//  Vector & F1,
	//  Vector & F2
	// ) const;

	bool
	defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const override { return true; }

	Distance
	atomic_interaction_cutoff() const override;

	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const override;

	bool
	use_extended_residue_pair_energy_interface() const override { return true; }

	void
	residue_pair_energy_ext(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		core::scoring::ResPairMinimizationData const & min_data,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;

private:

	void
	residue_pair_energy_one_way(
		conformation::Residue const & rsd1, // The other residue
		conformation::Residue const & rsd2, // The Mg(2+)
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;

	void
	eval_residue_pair(
		conformation::Residue const & ires,
		conformation::Residue const & jres,
		core::scoring::ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		core::scoring::EnergyMap & emap, // fill score values in here.
		core::scoring::EnergyMap const & weights, // for derivs.
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs) const;

	void
	eval_mg_interaction(
		conformation::Residue const & rsd1 /* other residue */,
		Size const atomno1 /*other atomno */,
		conformation::Residue const & rsd2 /* mg residue */,
		pose::Pose const & pose, // provides context
		core::scoring::EnergyMap & emap,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs /* other residue */,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs /* mg residue */
	) const;

	void
	eval_mg_residue_pair_derivatives(
		conformation::Residue const & rsd1 /* other residue */,
		conformation::Residue const & rsd2 /* mg residue */,
		pose::Pose const & pose, // provides context
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs /* other residue */,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs /* mg residue */
	) const;

	core::Size version() const override;

	core::scoring::magnesium::MgKnowledgeBasedPotentialOP mg_lig_knowledge_based_potential_;
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
	core::scoring::func::FuncCOP mg_sol_fade_func_;

};


} //scoring
} //core


#endif
