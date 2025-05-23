// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/P_AA_Energy.hh
/// @brief  Probability of observing an amino acid (NOT conditional on phi/psi), energy method declaration
/// @author Ron Jacak


#ifndef INCLUDED_core_energy_methods_P_AA_Energy_hh
#define INCLUDED_core_energy_methods_P_AA_Energy_hh

// Unit headers
#include <core/energy_methods/P_AA_Energy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/P_AA.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



class P_AA_Energy : public core::scoring::methods::ContextIndependentOneBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentOneBodyEnergy  parent;

public:

	P_AA_Energy();

	core::scoring::methods::EnergyMethodOP clone() const override;

	void residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const override;

	virtual
	Real
	eval_dof_derivative(
		id::DOF_ID const & dof_id, id::TorsionID const & tor_id, pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn, core::scoring::EnergyMap const & weights
	) const;

	/// @brief P_AA_Energy is context independent; indicates that no context graphs are required
	void indicate_required_context_graphs( utility::vector1< bool > & ) const override;


private:

	core::scoring::P_AA const & p_aa_;
	core::Size version() const override;

};

} // scoring
} // core


#endif // INCLUDED_core_energy_methods_P_AA_Energy_HH
