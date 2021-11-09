// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/FastDensEnergy.hh
/// @brief  Scoring a structure's fit to electron density
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_electron_density_FastDensEnergy_hh
#define INCLUDED_core_scoring_electron_density_FastDensEnergy_hh

#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers

namespace core {
namespace energy_methods {


class FastDensEnergy : public core::scoring::methods::ContextIndependentLRTwoBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentLRTwoBodyEnergy  parent;

public:

	/// constructor
	FastDensEnergy( core::scoring::methods::EnergyMethodOptions const & opts );


	/// clone
	core::scoring::methods::EnergyMethodOP clone() const override;

	/// lr container name
	core::scoring::methods::LongRangeEnergyType long_range_type() const override;


	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const override;

	bool
	defines_intrares_energy( core::scoring::EnergyMap const &  ) const override { return false; }

	/// scoring
	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap & emap
	) const override;

	void
	eval_intrares_energy(
		const core::conformation::Residue&,
		const core::pose::Pose&,
		const core::scoring::ScoreFunction&,
		core::scoring::EnergyMap&) const override {}

	/// derivatives
	void
	setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & sf) const override;

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

	///
	void
	finalize_total_energy(
		pose::Pose & ,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &
	) const override {}

	///  use the new minimizer interface
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override { return false; }

	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const override {};

	/// @details Has single atom energies, but not pairwise
	/// (Note that the residue-level calculation does this as a pairwise term for implementation reasons.
	/// The atomistic interface doesn't have the same limitations.)
	bool
	has_atomistic_energies() const override { return true; }

	void
	atomistic_energy(
		core::Size atmno, // Which atom?
		conformation::Residue const & rsd, // which residue?
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		scoring::EnergyMap & emap
	) const override;


private:
	bool scoreSymmComplex_;

	utility::vector1< core::Real > sc_scale_byres_;

	// helper function: make sure pose is setup for scoring
	bool
	pose_is_setup_for_density_scoring( pose::Pose const & pose) const;

	core::Size version() const override;
};


}
}

#endif

