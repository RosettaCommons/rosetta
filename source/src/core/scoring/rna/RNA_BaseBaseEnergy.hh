// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_BaseBaseEnergy.hh
/// @brief  Statistically derived base-base energy
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_BaseBaseEnergy_hh
#define INCLUDED_core_scoring_rna_RNA_BaseBaseEnergy_hh

// Unit Headers
#include <core/scoring/rna/RNA_BaseBaseEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers


namespace core {
namespace scoring {
namespace rna {

///

// Note that this could *almost* be a TwoBodyContextIndependentEnergy.
// The problem is that for low resolution structures,
// the base pair and base stack list needs to be filtered to
// disallow a single base edge from forming multiple base pairs, or
// for two bases to both pair and stack.
//


class RNA_BaseBaseEnergy : public methods::WholeStructureEnergy {
public:

	///
	RNA_BaseBaseEnergy();


	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap &// totals
	) const;

	virtual
	void
	eval_atom_derivative(
	  id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const &,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2 	) const;

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const {}


	// Undefinded, comented out to make python bindings complile
	//void
	//figure_out_rna_base_pairs_to_score() const;

	// Undefinded, comented out to make python bindings complile
	//void
	//figure_out_rna_base_stacks_to_score() const;

	// Undefinded, comented out to make python bindings complile
	//void
	//figure_out_which_interactions_to_score() const;


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	// const-ref to scoring database
	rna::RNA_LowResolutionPotential const & rna_low_resolution_potential_;

};


} //rna
} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
