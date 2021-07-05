// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNA_TorsionEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_TorsionEnergy_HH
#define INCLUDED_core_scoring_rna_RNA_TorsionEnergy_HH

// Unit Headers
#include <core/energy_methods/RNA_TorsionEnergy.fwd.hh>
#include <core/scoring/rna/RNA_EnergyMethodOptions.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/rna/RNA_TorsionPotential.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

namespace core {
namespace energy_methods {


class RNA_TorsionEnergy : public core::scoring::methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentTwoBodyEnergy  parent;
public:


	RNA_TorsionEnergy( core::scoring::rna::RNA_EnergyMethodOptions const & options,
		core::scoring::rna::RNA_TorsionPotentialOP rna_torsion_potential = nullptr );

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	/// @brief Evaluate the intra-residue constraint energy for a given residue
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
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

	bool
	defines_intrares_energy( core::scoring::EnergyMap const & ) const override { return true; }

	bool
	defines_residue_pair_energy( core::scoring::EnergyMap const & ) const { return true; }

	Distance
	atomic_interaction_cutoff() const override;

	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const override;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	core::scoring::rna::RNA_EnergyMethodOptions const & options_;
	core::scoring::rna::RNA_TorsionPotentialOP rna_torsion_potential_;

	core::Size version() const override;

};


} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
