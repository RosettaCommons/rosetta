// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_TorsionEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_TorsionEnergy_HH
#define INCLUDED_core_scoring_rna_RNA_TorsionEnergy_HH

// Unit Headers
#include <core/scoring/rna/RNA_TorsionEnergy.fwd.hh>
#include <core/scoring/rna/RNA_EnergyMethodOptions.fwd.hh>

// Package headers
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace rna {


class RNA_TorsionEnergy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;
public:


	RNA_TorsionEnergy( RNA_EnergyMethodOptions const & options,
										 RNA_TorsionPotentialOP rna_torsion_potential = 0 );

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	/// @brief Evaluate the intra-residue constraint energy for a given residue
	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const &,
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

	bool
	defines_intrares_energy( EnergyMap const & ) const { return true; }

	bool
	defines_residue_pair_energy( EnergyMap const & ) const { return true; }

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	RNA_EnergyMethodOptions const & options_;
	RNA_TorsionPotentialOP rna_torsion_potential_;

	virtual
	core::Size version() const;

};


} //rna
} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
