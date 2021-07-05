// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/GenericBondedEnergy.hh
/// @brief  A "generic" (atom-type-only-based) torsional potential
/// @author Hahnbeom Park and Frank DiMaio


#ifndef INCLUDED_core_energy_methods_GenericBondedEnergy_hh
#define INCLUDED_core_energy_methods_GenericBondedEnergy_hh

// Unit headers
#include <core/energy_methods/GenericBondedEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/GenericBondedPotential.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/MinimizationData.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#include <utility/vector1.hh>

namespace core {
namespace energy_methods {


class GenericBondedEnergy : public core::scoring::methods::ContextIndependentLRTwoBodyEnergy {
public:
	typedef core::scoring::methods::ContextIndependentLRTwoBodyEnergy  parent;

public:

	/// ctor
	GenericBondedEnergy( GenericBondedEnergy const & src );
	GenericBondedEnergy( core::scoring::methods::EnergyMethodOptions const & eopt );

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & sfxn ) const override;

	void
	setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & sfxn ) const override;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override { return false; }

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;

	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const override;

	bool
	defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const override { return true; }

	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &sfxn,
		core::scoring::EnergyMap & emap
	) const override;

	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		core::scoring::ResSingleMinimizationData const & /*res_data_cache*/,
		pose::Pose const & pose,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & atom_derivs
	) const override;

	bool
	defines_residue_pair_energy( pose::Pose const &, Size , Size ) const override { return ( true ); }

	bool
	defines_score_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		bool res_moving_wrt_eachother = false ) const override;

	core::scoring::methods::LongRangeEnergyType
	long_range_type() const override;

	/// @brief P_AA_pp_Energy is context independent; indicates that no
	/// context graphs are required
	void indicate_required_context_graphs( utility::vector1< bool > & ) const override;

	core::Size version() const override;

	virtual
	Distance
	atomic_interaction_cutoff() const { return 0.0; }

	bool score_full() const { return score_full_; }
	bool score_hybrid() const { return score_hybrid_; }

private:
	core::scoring::GenericBondedPotential const &potential_;
	bool score_full_;
	bool score_hybrid_;
};

} // scoring
} // core


#endif // INCLUDED_core_scoring_GenericBondedEnergy_hh

