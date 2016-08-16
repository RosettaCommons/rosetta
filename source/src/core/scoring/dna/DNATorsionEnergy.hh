// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/dna/DNATorsionEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Rhiju Das
/// @author Jim Havranek


#ifndef INCLUDED_core_scoring_dna_DNATorsionEnergy_HH
#define INCLUDED_core_scoring_dna_DNATorsionEnergy_HH

// Unit Headers
#include <core/scoring/dna/DNATorsionEnergy.fwd.hh>

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
namespace scoring {
namespace dna {


class DNATorsionEnergy : public methods::ContextIndependentTwoBodyEnergy {
public:
	typedef ContextIndependentTwoBodyEnergy parent;
public:


	DNATorsionEnergy();

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	setup_for_scoring( pose::Pose &pose, ScoreFunction const &scfxn ) const;

	// call the cst setup_for_derivatives wrapper
	virtual
	void
	setup_for_derivatives( pose::Pose &pose, ScoreFunction const &scfxn ) const;


	// virtual
	// void
	// setup_for_packing( pose::Pose & pose, pack::task::PackerTask const & ) const;

	/// @brief Evaluate the intra-residue constraint energy for a given residue
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
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	/// called at the end of energy evaluation
	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

	/// called at the end of energy evaluation
	virtual
	void
	finalize_after_derivatives(
		pose::Pose & pose,
		ScoreFunction const &
	) const;


	/// called during gradient-based minimization inside dfunc
	/**
	F1 and F2 are not zeroed -- contributions from this atom are
	just summed in
	**/
	virtual
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &, // domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;


	/// uses the dof constraints
	Real
	eval_dof_derivative(
		id::DOF_ID const & id,
		id::TorsionID const & tor,
		pose::Pose const & pose,
		ScoreFunction const & scorefxn,
		EnergyMap const & weights
	) const;


	bool
	defines_intrares_energy( EnergyMap const & ) const { return true; }

	bool
	defines_residue_pair_energy( EnergyMap const & ) const { return true; }

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const;

	virtual
	core::Size version() const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	dna::DNATorsionPotential const & dna_torsion_potential_;

	mutable constraints::ConstraintSetOP dna_torsion_constraints_;
	mutable constraints::ConstraintSetOP dna_sugar_close_constraints_;
	mutable constraints::ConstraintSetOP dna_base_distance_constraints_;

	mutable bool constraints_ready_;
	bool verbose_;

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
