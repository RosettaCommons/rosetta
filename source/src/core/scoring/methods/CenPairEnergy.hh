// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/CenPairEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_methods_CenPairEnergy_hh
#define INCLUDED_core_scoring_methods_CenPairEnergy_hh

// Unit Headers

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/EnvPairPotential.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class CenPairEnergy : public ContextIndependentTwoBodyEnergy  {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;
public:


	CenPairEnergy();


	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const &, // pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;


	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap &// totals
	) const;

	/// This method *should* admit to defining intraresidue energies
	virtual
	bool
	defines_intrares_energy( EnergyMap const & ) const { return false; }

	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {}

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const {}


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	// const-ref to scoring database
	EnvPairPotential const & potential_;
	virtual
	core::Size version() const;

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
