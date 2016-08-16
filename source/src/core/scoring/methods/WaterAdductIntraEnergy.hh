// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/WaterAdductIntraEnergy.hh
/// @brief  Energetic offset/cost for placing a water adduct on an amino or nucleic acid
/// @author Jim Havranek


#ifndef INCLUDED_core_scoring_methods_WaterAdductIntraEnergy_hh
#define INCLUDED_core_scoring_methods_WaterAdductIntraEnergy_hh

// Unit headers
#include <core/scoring/methods/WaterAdductIntraEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class WaterAdductIntraEnergy : public ContextIndependentOneBodyEnergy  {
public:
	typedef ContextIndependentOneBodyEnergy  parent;

public:

	/// ctor
	WaterAdductIntraEnergy();

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


	virtual
	Real
	eval_dof_derivative(
		id::DOF_ID const & dof_id,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights
	) const;

	/// @brief WaterAdductIntraEnergy is context independent; indicates that no
	/// context graphs are required
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;
	virtual
	core::Size version() const;

	// data
private:

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_methods_WaterAdductIntraEnergy_HH
