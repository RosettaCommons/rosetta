// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNA_VDW_Energy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_rna_RNA_VDW_Energy_hh
#define INCLUDED_core_scoring_rna_RNA_VDW_Energy_hh

// Unit Headers
#include <core/scoring/rna/RNA_VDW_Energy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/rna/RNA_AtomVDW.fwd.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace scoring {
namespace rna {


class RNA_VDW_Energy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;

public:


	RNA_VDW_Energy();

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
	setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;

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
	) const {}


	virtual
	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const & scorefxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

private:

	Size
	get_vdw_atom_number(
		utility::vector1< utility::vector1< Size > > const & atom_numbers_for_vdw_calculation,
		Size const & pos1,
		Size const & i ) const;

	Size
	get_vdw_atom_number(
		char const which_nucleotide,
		Size const & i ) const;

	Real
	evaluate_rnp_vdw_score( 
		conformation::Residue const & rna_rsd,
		conformation::Residue const & protein_rsd,
		rna::RNA_ScoringInfo const & rna_scoring_info,
		bool const & centroid_mode ) const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	//Helper functions just to get things set up.
	void
	setup_atom_numbers_for_vdw_calculation( pose::Pose & pose ) const;

	// Data!
	// const-ref to scoring database
	rna::RNA_AtomVDW const & rna_atom_vdw_;

	Real const vdw_scale_factor_;
	virtual
	core::Size version() const;
};


} //rna
} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
