// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_JR_SuiteEnergy.hh
/// @brief  Jane Richardson style RNA Suite Energy
/// @author Fang-Chieh Chou

#ifndef INCLUDED_core_scoring_rna_RNA_JR_SuiteEnergy_hh
#define INCLUDED_core_scoring_rna_RNA_JR_SuiteEnergy_hh

// Unit Headers
#include <core/scoring/rna/RNA_JR_SuiteEnergy.fwd.hh>

// Package headers
#include <core/conformation/Atom.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/rna/RNA_SuiteName.hh>

// Utility headers
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>

namespace core {
namespace scoring {
namespace rna {

///
class RNA_JR_SuiteEnergy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy parent;

	RNA_JR_SuiteEnergy();

	/// clone
	virtual	methods::EnergyMethodOP	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	/// @brief Evaluate the intra-residue constraint energy for a given residue
	virtual	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {};

	///
	virtual	void
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
	virtual	void
	eval_atom_derivative(
		id::AtomID const &,
		pose::Pose const &,
		kinematics::DomainMap const &,
		ScoreFunction const &,
		EnergyMap const &,
		Vector &,
		Vector &
	) const {};

	bool defines_intrares_energy( EnergyMap const & ) const { return true; }

	bool
	defines_residue_pair_energy( EnergyMap const & ) const { return true; }

	virtual	Distance atomic_interaction_cutoff() const { return 0; }

	virtual	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const {};
/////////////////////////////////////////////////////////////////////////////
// data
/////////////////////////////////////////////////////////////////////////////

private:
	virtual	Size version() const { return 0; }

	core::pose::rna::RNA_SuiteName suitename_;
};

} //rna
} //scoring
} //core

#endif // INCLUDED_core_scoring_methods_RNA_LJ_BaseEnergy_HH
