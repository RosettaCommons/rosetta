// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ElecDensEnergy.hh
/// @brief  Scoring a structure's fit to electron density
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_electron_density_ElecDensCenEnergy_hh
#define INCLUDED_core_scoring_electron_density_ElecDensCenEnergy_hh


// Unit Headers

// Package headers
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers

namespace core {
namespace energy_methods {


class ElecDensCenEnergy : public core::scoring::methods::ContextDependentLRTwoBodyEnergy  {
public:
	typedef core::scoring::methods::ContextDependentLRTwoBodyEnergy  parent;

public:


	ElecDensCenEnergy();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	core::scoring::methods::LongRangeEnergyType
	long_range_type() const override;


	bool defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const override;

	void setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	void setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & sf) const override;

	bool defines_intrares_energy( core::scoring::EnergyMap const &  ) const override { return true; }

	/// @brief Evaluate the intra-residue constraint energy for a given residue
	void eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap & emap
	) const override ;


	void residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap & emap
	) const override;

	using core::scoring::methods::ContextDependentLRTwoBodyEnergy::finalize_total_energy;

	/// called at the end of energy evaluation
	virtual void finalize_total_energy(
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const;


	/// called during gradient-based minimization inside dfunc
	void eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &, // domain_map,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;


	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const override {};

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////


private:
	mutable bool pose_is_proper;
	mutable core::Real structure_score;
	mutable int nreses;
	core::Size version() const override;
};


}
}

#endif

