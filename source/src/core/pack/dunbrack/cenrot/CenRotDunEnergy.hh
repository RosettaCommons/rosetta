// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/CenRotDunEnergy.cc
/// @brief  CenRot version of centroid Dunbrack energy
/// @author Yuan Liu


#ifndef INCLUDED_core_pack_dunbrack_CenRotDunEnergy_hh
#define INCLUDED_core_pack_dunbrack_CenRotDunEnergy_hh

// Unit Headers

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace pack {
namespace dunbrack {
namespace cenrot {


class CenRotDunEnergy : public scoring::methods::ContextIndependentOneBodyEnergy  {
public:
	typedef scoring::methods::ContextIndependentOneBodyEnergy  parent;
public:


	CenRotDunEnergy();


	/// clone
	virtual
	scoring::methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	setup_for_scoring( pose::Pose & pose, scoring::ScoreFunction const & ) const;

	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		scoring::EnergyMap & emap
	) const;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	virtual
	void
	eval_residue_derivatives(
		conformation::Residue const & rsd,
		scoring::ResSingleMinimizationData const & min_data,
		pose::Pose const & pose,
		scoring::EnergyMap const & weights,
		utility::vector1< scoring::DerivVectorPair > & atom_derivs
	) const;

	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		scoring::ScoreFunction const &,
		scoring::EnergyMap &// totals
	) const;

	/// @brief Yes.  The DunbrackEnergy defines derivatives
	/// for phi/psi and the chi dihedrals.
	virtual
	bool
	defines_dof_derivatives( pose::Pose const & p ) const;

	/// @brief Evaluate the phi/psi and chi dihedral derivatives
	/// for the input residue.
	Real
	eval_residue_dof_derivative(
		conformation::Residue const & rsd,
		scoring::ResSingleMinimizationData const & min_data,
		id::DOF_ID const & dof_id,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		scoring::EnergyMap const & weights
	) const;

	/// @brief Deprecated.
	virtual
	Real
	eval_dof_derivative(
		id::DOF_ID const & dof_id,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		scoring::EnergyMap const & weights
	) const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const {}


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	virtual
	core::Size version() const;

};


}
}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
