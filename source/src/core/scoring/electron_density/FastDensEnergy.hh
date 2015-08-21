// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
namespace scoring {
namespace electron_density {


class FastDensEnergy : public methods::ContextIndependentLRTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentLRTwoBodyEnergy  parent;

public:

	/// constructor
	FastDensEnergy( methods::EnergyMethodOptions const & opts );


	/// clone
	virtual methods::EnergyMethodOP clone() const;

	/// lr container name
	methods::LongRangeEnergyType long_range_type() const;


	virtual bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const &  ) const { return false; }

	/// scoring
	virtual void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	virtual void
	eval_intrares_energy(
		const core::conformation::Residue&,
		const core::pose::Pose&,
		const core::scoring::ScoreFunction&,
		core::scoring::EnergyMap&) const {}

	/// derivatives
	virtual void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sf) const;

	virtual void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const & min_data,
		pose::Pose const &,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;


	virtual void
	finalize_total_energy(
		pose::Pose & ,
		ScoreFunction const &,
		EnergyMap &
	) const {}

	///  use the new minimizer interface
	virtual bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const {};

private:
	bool scoreSymmComplex_;

	utility::vector1< core::Real > sc_scale_byres_;

	// helper function: make sure pose is setup for scoring
	bool
	pose_is_setup_for_density_scoring( pose::Pose const & pose) const;

	virtual
	core::Size version() const;
};


}
}
}

#endif

