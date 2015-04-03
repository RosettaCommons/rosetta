// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
	virtual ~DunbrackEnergy();

	/// clone
	virtual
	scoring::methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		scoring::EnergyMap & emap
	) const;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

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

	/// @brief DunbrackEnergy is context independent; indicates that no
	/// context graphs are required
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;


	// data
private:
	// changed to & so that it does not duplicate the dunbrack lib; arguably should be COP.
	//RotamerLibrary const & rot_lib_;

	//mutable dunbrack::RotamerLibraryScratchSpaceOP scratch_;
virtual
core::Size version() const;

};

} // dunbrack
} // pack
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
