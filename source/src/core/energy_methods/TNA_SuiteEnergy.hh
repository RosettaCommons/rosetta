// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/TNA_SuiteEnergy.hh
/// @brief  TNA Suite Energy
/// @author Andy Watkins

#ifndef INCLUDED_core_scoring_rna_TNA_SuiteEnergy_hh
#define INCLUDED_core_scoring_rna_TNA_SuiteEnergy_hh

// Unit Headers
#include <core/energy_methods/TNA_SuiteEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/rna/TNA_SuitePotential.fwd.hh>
#include <core/scoring/rna/RNA_EnergyMethodOptions.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

namespace core {
namespace energy_methods {

class TNA_SuiteEnergy : public core::scoring::methods::ContextIndependentLRTwoBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentLRTwoBodyEnergy parent;

	TNA_SuiteEnergy( core::scoring::rna::RNA_EnergyMethodOptions const & options );

	// clone
	core::scoring::methods::EnergyMethodOP clone() const override { return utility::pointer::make_shared< TNA_SuiteEnergy >( options_ ); }

	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	// scoring
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
	) const override {}

	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResPairMinimizationData const & min_data,
		pose::Pose const &,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const override;

	bool
	defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const override { return false; }

	bool
	defines_intrares_dof_derivatives( pose::Pose const & ) const override { return false; }

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override { return false; }

	bool
	defines_residue_pair_energy(
		pose::Pose const &,
		Size res1,
		Size res2 ) const override
	{
		return ( res1 == (res2+1) || res2 == (res1+1) );
	}

	virtual Distance atomic_interaction_cutoff() const { return 0; }

	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const override {};

	core::scoring::methods::LongRangeEnergyType
	long_range_type() const override { return core::scoring::methods::tna_suite_lr; }

private:
	Size version() const override { return 420; }

	bool get_f1_f2(
		id::TorsionID const & torsion_id,
		pose::Pose const & pose,
		utility::vector1< id::AtomID > & atom_ids,
		utility::vector1< Vector > & f1s,
		utility::vector1< Vector > & f2s
	) const;

	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		Real const & weights,
		utility::vector1< core::scoring::DerivVectorPair> & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair> & r2_atom_derivs,
		core::scoring::rna::TNA_SuitePotentialCOP rna_suite_potential ) const;

private:

	core::scoring::rna::RNA_EnergyMethodOptions const & options_;
	core::scoring::rna::TNA_SuitePotentialOP tna_suite_potential_;

};

} //scoring
} //core

#endif // INCLUDED_core_energy_methods_TNA_SuiteEnergy_HH
