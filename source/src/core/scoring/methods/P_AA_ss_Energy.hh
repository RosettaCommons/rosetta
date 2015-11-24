// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/P_AA_ss_Energy.hh
/// @brief  Probability of observing an amino acid (NOT conditional on phi/psi), energy method declaration
/// @author Ron Jacak


#ifndef INCLUDED_core_scoring_methods_P_AA_ss_Energy_hh
#define INCLUDED_core_scoring_methods_P_AA_ss_Energy_hh

// Unit headers
#include <core/scoring/methods/P_AA_ss_Energy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/P_AA_ss.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

// fd: this probably should be context-dependent
class P_AA_ss_Energy : public ContextIndependentOneBodyEnergy  {
public:
	typedef ContextIndependentOneBodyEnergy  parent;

public:

	P_AA_ss_Energy();

	virtual
	EnergyMethodOP clone() const;

	virtual
	void
	setup_for_scoring(
		pose::Pose & pose,
		ScoreFunction const & ) const;

	virtual
	void
	setup_for_minimizing(
		pose::Pose & ,
		ScoreFunction const & ,
		kinematics::MinimizerMapBase const &
	) const;

	virtual
	Real
	eval_dof_derivative(
		id::DOF_ID const & dof_id, id::TorsionID const & tor_id, pose::Pose const & pose,
		ScoreFunction const & sfxn, EnergyMap const & weights
	) const;

	virtual
	void residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief P_AA_ss_Energy is context independent; indicates that no context graphs are required
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;

private:

	P_AA_ss const & P_AA_ss_;

	virtual
	core::Size version() const;
};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_methods_P_AA_ss_Energy_HH
