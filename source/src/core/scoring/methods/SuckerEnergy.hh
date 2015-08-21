// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SuckerEnergy.hh
/// @brief  energy attracting atoms to sucker atoms
/// @author Will Sheffler


#ifndef INCLUDED_core_scoring_methods_SuckerEnergy_hh
#define INCLUDED_core_scoring_methods_SuckerEnergy_hh

// Unit Headers
#include <core/scoring/methods/SuckerEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <numeric/interpolation/spline/Interpolator.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace scoring {
namespace methods {


class SuckerEnergy : public ContextIndependentTwoBodyEnergy  {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;
public:


	SuckerEnergy();
	~SuckerEnergy();


	/// clone
	virtual
	EnergyMethodOP
	clone() const;


	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & scorefxn,
		EnergyMap & emap
	) const;


	virtual
	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;


	virtual
	Distance
	atomic_interaction_cutoff() const;

	/// @details non-virtual accessor for speed
	Distance
	interaction_cutoff() const;


	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const {
		return true;
	}

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	numeric::interpolation::spline::InterpolatorOP interp_;
	virtual
	core::Size version() const;

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
