// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/DNATorsionEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Rhiju Das
/// @author Jim Havranek


#ifndef INCLUDED_core_energy_methods_DNATorsionEnergy_HH
#define INCLUDED_core_energy_methods_DNATorsionEnergy_HH

// Unit Headers
#include <core/energy_methods/DNATorsionEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
//#include <core/scoring/dna/DNATorsionPotential.fwd.hh>
#include <core/scoring/dna/DNATorsionPotential.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
//#include <core/pack/task/PackerTask.hh>

#include <map>

// Utility headers


namespace core {
namespace energy_methods {


class DNATorsionEnergy : public core::scoring::methods::ContextIndependentTwoBodyEnergy {
public:
	typedef core::scoring::methods::ContextIndependentTwoBodyEnergy parent;

	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::scoring::ScoreFunction ScoreFunction;
public:


	DNATorsionEnergy();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	setup_for_scoring( pose::Pose &pose, core::scoring::ScoreFunction const &scfxn ) const override;

	// call the cst setup_for_derivatives wrapper
	void
	setup_for_derivatives( pose::Pose &pose, core::scoring::ScoreFunction const &scfxn ) const override;


	// virtual
	// void
	// setup_for_packing( pose::Pose & pose, pack::task::PackerTask const & ) const;

	/// @brief Evaluate the intra-residue constraint energy for a given residue
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap & emap
	) const override;


	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;

	/// called at the end of energy evaluation
	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const override;

	/// called at the end of energy evaluation
	void
	finalize_after_derivatives(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &
	) const override;


	/// called during gradient-based minimization inside dfunc
	/**
	F1 and F2 are not zeroed -- contributions from this atom are
	just summed in
	**/
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &, // domain_map,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;


	/// uses the dof constraints
	Real
	eval_dof_derivative(
		id::DOF_ID const & id,
		id::TorsionID const & tor,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & scorefxn,
		core::scoring::EnergyMap const & weights
	) const;


	bool
	defines_intrares_energy( core::scoring::EnergyMap const & ) const override { return true; }

	bool
	defines_residue_pair_energy( core::scoring::EnergyMap const & ) const { return true; }

	Distance
	atomic_interaction_cutoff() const override;

	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const override;

	core::Size version() const override;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	core::scoring::dna::DNATorsionPotential const & dna_torsion_potential_;

	mutable core::scoring::constraints::ConstraintSetOP dna_torsion_constraints_;
	mutable core::scoring::constraints::ConstraintSetOP dna_sugar_close_constraints_;
	mutable core::scoring::constraints::ConstraintSetOP dna_base_distance_constraints_;

	mutable bool constraints_ready_;
	bool verbose_;

};

}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
