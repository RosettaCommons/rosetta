// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/methods/carbohydrates/SugarBackboneEnergy.hh
/// @brief   Method declarations and simple accessor definitions for SugarBackboneEnergy.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_scoring_methods_carbohydrates_SugarBackboneEnergy_HH
#define INCLUDED_core_scoring_methods_carbohydrates_SugarBackboneEnergy_HH

// Unit header
#include <core/scoring/methods/carbohydrates/SugarBackboneEnergy.fwd.hh>
#include <core/scoring/carbohydrates/CHIEnergyFunction.fwd.hh>

// Package headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>

// Project headers
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <map>

namespace core {
namespace scoring {
namespace methods {
namespace carbohydrates {


/// @details  This class makes use of the "CarboHydrate Intrinsic" (CHI) energy function developed by Woods Lab.\n
/// Carbohydrate phi angles are scored based on whether they are at alpha or beta linkages.\n
/// Psi angles are scored based on whether they are at ->2-axial, ->3-equatorial, or ->4-axial OR ->2-equatorial,
/// ->3-axial, or ->4-equatorial linkages.\n
/// ->6 linkages (with omega angles) are not scored.
/// @ref      A.K. Nivedha et al. J. Comput. Chem. 2014, 35, 526-39
class SugarBackboneEnergy : public ContextIndependentOneBodyEnergy {
public:  // Standard Methods //////////////////////////////////////////////////
	/// @brief  Default constructor
	SugarBackboneEnergy();


public:  // General EnergyMethod Methods //////////////////////////////////////
	virtual EnergyMethodOP clone() const;

	/// @brief   Should this EnergyMethod have score and derivative evaluation evaluated ONLY in the context of a whole
	/// Pose?
	/// @return  false
	virtual bool minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	/// @brief    Indicate in the context-graphs-required list which context-graphs this energy method requires that the
	/// Pose maintains when doing neighbor evaluation.
	/// @details  not implemented for SugarBackboneEnergy
	virtual void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const {}


public:  // OneBodyEnergy Methods /////////////////////////////////////////////
	/// @brief  Evaluate the one-body carbohydrate backbone energies for a particular residue, in the context of the
	/// given Pose, and increment those energies in the input Emap.
	virtual void residue_energy( conformation::Residue const & rsd, pose::Pose const & pose, EnergyMap & emap ) const;

	/// @brief   Should the dof_derivative interface be used for this EnergyMethod when calculating derivatives?
	/// @return  true
	virtual bool defines_dof_derivatives( pose::Pose const & /* pose */ ) const { return true; }

	/// @brief    Evaluate the DoF derivative for a particular residue.
	virtual core::Real eval_residue_dof_derivative(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		id::DOF_ID const & dof_id,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose,
		ScoreFunction const & sf,
		EnergyMap const & weights ) const;


private:  // Private methods //////////////////////////////////////////////////
	virtual core::Size version() const { return 1; }  // initial versioning

private:  // Private Data /////////////////////////////////////////////////////
	// the "CarboHydrate Intrinsic" (CHI) energy function developed by Woods Lab
	scoring::carbohydrates::CHIEnergyFunction const & E_;



};

}  // namespace carbohydrates
}  // namespace methods
}  // namespace scoring
}  // namespace core

#endif  // INCLUDED_core_scoring_methods_carbohydrates_SugarBackboneEnergy_HH
