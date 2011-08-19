// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/.cc
///
/// @brief
/// @author Ian W. Davis
/// @author Kristian Kaufmann


#include <core/scoring/methods/CSD_TorsionEnergy.hh>

#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <numeric/xyz.functions.hh>

namespace core {
namespace scoring {
namespace methods {

/// @brief THIS CLASS REQUIRES AN EnergyMethodCreator TO WORK PROPERLY
CSD_TorsionEnergy::CSD_TorsionEnergy():
	ContextIndependentTwoBodyEnergy()
{}


CSD_TorsionEnergy::~CSD_TorsionEnergy() {}


EnergyMethodOP
CSD_TorsionEnergy::clone() const
{
	return new CSD_TorsionEnergy();
}


void
CSD_TorsionEnergy::setup_for_packing( pose::Pose & pose, pack::task::PackerTask const & ) const
{
	pose.update_residue_neighbors();
}


void
CSD_TorsionEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


void
CSD_TorsionEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


bool
CSD_TorsionEnergy::defines_intrares_energy( EnergyMap const & ) const
{
	return true;
}


void
CSD_TorsionEnergy::residue_pair_energy(
	 conformation::Residue const & ,
	 conformation::Residue const & ,
	 pose::Pose const & ,
	 ScoreFunction const & ,
	 EnergyMap & emap
) const
{
	// add energy to emap (unweighted)
	emap[ csd_torsion ] += 0;
}


void
CSD_TorsionEnergy::eval_intrares_energy(
	conformation::Residue const & ,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap &
) const
{
}


Real
CSD_TorsionEnergy::eval_dof_derivative(
	id::DOF_ID const &,// dof_id,
	id::TorsionID const & ,
	pose::Pose const & ,
	ScoreFunction const &,// sfxn,
	EnergyMap const & weights
) const
{
	// note that the atomtree Oomega dofs are in radians
	// use degrees since dE/dangle has angle in denominator
	Real deriv = 0;
	return numeric::conversions::degrees( weights[ csd_torsion ] * deriv );
}


Distance
CSD_TorsionEnergy::atomic_interaction_cutoff() const
{
	return 0.0;
}


void
CSD_TorsionEnergy::indicate_required_context_graphs(utility::vector1< bool > & ) const
{}


} // namespace methods
} // namespace scoring
} // namespace core
