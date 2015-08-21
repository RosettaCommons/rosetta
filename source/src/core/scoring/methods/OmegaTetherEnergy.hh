// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/OmegaTetherEnergy.hh
/// @brief  OmegaTether energy method class declaration
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_methods_OmegaTetherEnergy_hh
#define INCLUDED_core_scoring_methods_OmegaTetherEnergy_hh

// Unit headers
#include <core/scoring/methods/OmegaTetherEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/OmegaTether.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class OmegaTetherEnergy : public ContextIndependentOneBodyEnergy  {
public:
	typedef ContextIndependentOneBodyEnergy  parent;
public:

	/// ctor
	OmegaTetherEnergy();

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

	/// @brief Use the dof_derivative interface for this energy method when
	/// calculating derivatives?  It is possible to define both dof_derivatives and
	/// atom-derivatives; they are not mutually exclusive.
	virtual
	bool
	defines_dof_derivatives( pose::Pose const & p ) const;

	/// @brief Evaluate the DOF derivative for a particular residue.  The Pose merely serves as context,
	/// and the input residue is not required to be a member of the Pose.
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


	virtual
	Real
	old_eval_dof_derivative(
		id::DOF_ID const &,// dof_id,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		ScoreFunction const &,// sfxn,
		EnergyMap const & weights
	) const;

	/// @brief OmegaTether Energy is context independent and thus indicates that no context graphs need to
	/// be maintained by class Energies
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const;

	// data
private:
	OmegaTether const & potential_;
	virtual
	core::Size version() const;

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
