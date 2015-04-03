// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/scoring/methods/SpecialRotamerEnergy.hh
/// @brief  Adds a bonus to any rotamer that is flagged
/// @author sthyme, sthyme@gmail.com, Feb 2010


#ifndef INCLUDED_protocols_scoring_methods_SpecialRotamerEnergy_hh
#define INCLUDED_protocols_scoring_methods_SpecialRotamerEnergy_hh

// Unit headers
#include <protocols/scoring/methods/SpecialRotamerEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace scoring {
namespace methods {


class SpecialRotamerEnergy : public core::scoring::methods::ContextIndependentOneBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentOneBodyEnergy  parent;

public:

	/// ctor
	SpecialRotamerEnergy();

	/// clone
	virtual
	core::scoring::methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	virtual
	void
	residue_energy(
		core::conformation::Residue const & rsd,
		core::pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;


	virtual
	core::Real
	eval_dof_derivative(
		core::id::DOF_ID const & dof_id,
		core::id::TorsionID const & tor_id,
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights
	) const;

	/// @brief SpecialRotamerEnergy is context independent; indicates that no
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
} // protocols


#endif // INCLUDED_protocols_scoring_methods_SpecialRotamerEnergy_HH
