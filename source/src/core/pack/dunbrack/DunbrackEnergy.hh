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

private:
	// methods

	/// @brief Given a mainchain torsion index and a ResidueType, get the index of the corresponding torsion in the
	/// data stored in the Dunbrack scratch space.
	/// @details For most residue types, this just returns torsion_index.  The index is only different in cases in which
	/// a residue type has rotamers that depend on a subset of mainchain torsions.  For example, if a residue's rotamers
	/// depended on mainchain torsions 2, 3, and 4, then the scratch indices 1, 2, and 3 would correspond to mainchain
	/// torsions 2, 3, and 4, respectively.  This function returns 0 if torsion_index is a torsion on which rotamers do
	/// not depend.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	core::Size get_scratch_index( core::Size const torsion_index, core::conformation::Residue const &rsd ) const;

	virtual
	core::Size version() const;

	// data
private:


};

} // dunbrack
} // pack
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
