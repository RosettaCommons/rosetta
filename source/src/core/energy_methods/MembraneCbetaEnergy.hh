// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/MembraneCbetaEnergy.hh
/// @brief  Membrane Environment Cbeta Energy
/// @author Bjorn Wallner


#ifndef INCLUDED_core_scoring_methods_MembraneCbetaEnergy_hh
#define INCLUDED_core_scoring_methods_MembraneCbetaEnergy_hh

// Unit Headers
//#include <core/energy_methods/MembraneCbetaEnergy.fwd.hh>
// you cannot #include yourself #include <core/energy_methods/MembraneCbetaEnergy.hh>
#include <core/scoring/MembraneTopology.fwd.hh>
// Package headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/MembranePotential.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace scoring {
namespace methods {


class MembraneCbetaEnergy : public ContextDependentOneBodyEnergy  {
public:
	typedef ContextDependentOneBodyEnergy  parent;

public:


	MembraneCbetaEnergy();


	/// clone
	EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const override;

	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const override;


	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap &// totals
	) const override;

	void indicate_required_context_graphs( utility::vector1< bool > & ) const override {}

	//bool allow_scoring(Size const seqpos) const;
	MembraneTopology const & MembraneTopology_from_pose( pose::Pose const & pose ) const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	// const-ref to scoring database
	MembranePotential const & potential_;
	core::Size version() const override;

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
