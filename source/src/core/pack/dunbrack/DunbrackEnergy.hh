// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/dunbrack/DunbrackEnergy.hh
/// @brief
/// @author


#ifndef INCLUDED_core_pack_dunbrack_DunbrackEnergy_hh
#define INCLUDED_core_pack_dunbrack_DunbrackEnergy_hh

// Unit headers
#include <core/pack/dunbrack/DunbrackEnergy.fwd.hh>

// Package headers
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace dunbrack {


class DunbrackEnergy : public scoring::methods::ContextIndependentOneBodyEnergy  {
public:
	typedef ContextIndependentOneBodyEnergy  parent;
	typedef dunbrack::RotamerLibrary RotamerLibrary;

public:

	/// @brief ctor
	DunbrackEnergy();

	/// @brief dtor
	~DunbrackEnergy() override;

	/// clone
	scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		scoring::EnergyMap & emap
	) const override;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override { return false; }

	/// @brief Yes.  The DunbrackEnergy defines derivatives
	/// for phi/psi and the chi dihedrals.
	bool
	defines_dof_derivatives( pose::Pose const & p ) const override;

	/// @brief Return a vector of the atoms that determine the DOFs covered by
	/// the rotamer libraries; these may extend beyond the residue in question
	/// and therefore must be described using PartialAtomIDs.
	utility::vector1< id::PartialAtomID >
	atoms_with_dof_derivatives( conformation::Residue const & res, pose::Pose const & ) const override;

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
	) const override;

	/// @brief DunbrackEnergy is context independent; indicates that no
	/// context graphs are required
	void indicate_required_context_graphs( utility::vector1< bool > & ) const override;

private:
	// methods


	core::Size version() const override;

	// data
private:


};

} // dunbrack
} // pack
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
