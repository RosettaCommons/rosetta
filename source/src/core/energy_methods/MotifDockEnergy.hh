// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/MotifDockEnergy.hh
/// @brief  Adaptation of Motif score for Docking
/// @author Nick Marze (nickmarze@gmail.com)


#ifndef INCLUDED_core_scoring_methods_MotifDockEnergy_hh
#define INCLUDED_core_scoring_methods_MotifDockEnergy_hh

// Unit Headers

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/energy_methods/MotifDockEnergy.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <numeric/xyzTransform.hh>

#include <utility/vector1.hh>


typedef numeric::xyzTransform<core::Real> Xform;


namespace core {
namespace scoring {
namespace methods {


class MotifDockEnergy : public ContextIndependentTwoBodyEnergy  {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;
public:


	MotifDockEnergy();


	/// clone
	EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const override;

	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const override;

	bool
	defines_intrares_energy( EnergyMap const & weights ) const override;

	void indicate_required_context_graphs( utility::vector1< bool > & ) const override {};

	core::Size version() const override;

	Distance
	atomic_interaction_cutoff() const override;
};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
