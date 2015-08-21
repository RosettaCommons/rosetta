// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RamachandranEnergy.hh
/// @brief  Ramachandran energy method class declaration
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_methods_RamachandranEnergy_hh
#define INCLUDED_core_scoring_methods_RamachandranEnergy_hh

// Unit headers
#include <core/scoring/methods/RamachandranEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/Ramachandran.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class RamachandranEnergy : public ContextIndependentOneBodyEnergy  {
public:
	typedef ContextIndependentOneBodyEnergy  parent;
public:

	/// ctor
	RamachandranEnergy();

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	/// @brief The ramachandran energy defines derivatives for protein backbone torsion angles
	virtual
	bool
	defines_dof_derivatives( pose::Pose const & p ) const;

	/// @brief Evaluate the phi or psi derivative for a particular residue
	virtual
	Real
	eval_residue_dof_derivative(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		id::DOF_ID const & dof_id,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights
	) const;

	/// @brief NOTE: non-virtual function interface.
	Real
	eval_dof_derivative(
		id::DOF_ID const &,// dof_id,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		ScoreFunction const &,// sfxn,
		EnergyMap const & weights
	) const;

	/// @brief Ramachandran Energy is context independent and thus indicates that no context graphs need to
	/// be maintained by class Energies
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const;

	// data
private:
	Ramachandran const & potential_;
	virtual
	core::Size version() const;

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
