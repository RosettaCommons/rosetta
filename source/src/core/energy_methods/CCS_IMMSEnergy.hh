// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/energy_methods/CCS_IMMSEnergy.hh
/// @brief  energy term use for scoring predicted CCS_IMMSEnergy
/// @brief  this energy method does not return derivatives, and is not compatible with the minimizer
/// @brief  core::Real parcs_ccs (.., .., ..) needs to be non-member function because CCS_IMMSEnergy() requires Experimental CCS, but ParcsCCS calculation does not. So defining it outside the CCS_IMMSEnergy
/// @author SM Baargeen Alam Turzo <turzo.1@osu.edu>

#ifndef INCLUDED_core_energy_methods_CCS_IMMSEnergy_hh
#define INCLUDED_core_energy_methods_CCS_IMMSEnergy_hh

#include <core/energy_methods/CCS_IMMSEnergy.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

#include <core/energy_methods/CCS_IMMSEnergyCreator.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>



#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace energy_methods {


/////////////////////////////////////////////////////////////////////////////
// CCS Calculation
/////////////////////////////////////////////////////////////////////////////

core::Real
parcs_ccs(core::pose::Pose &mypose, core::Size const nrot, core::Real const prad);

class CCS_IMMSEnergy : public core::scoring::methods::WholeStructureEnergy {
public:
	typedef core::scoring::methods::WholeStructureEnergy parent;

public:
	CCS_IMMSEnergy();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	core::Real
	calc_IMMS_score(const core::Real &CCS_pred, const core::Real &CCS_exp) const;


	void
	finalize_total_energy(core::pose::Pose &mypose, core::scoring::ScoreFunction const &, core::scoring::EnergyMap & emap) const override;


	void indicate_required_context_graphs( utility::vector1< bool > & ) const override;
	core::Size version() const override;


private:
	core::Real ccs_exp_;
	core::Size nrot_=300;
	core::Real prad_=1.0;

};

} // scoring
} // core

#endif // INCLUDED_core_energy_methods_CCS_IMMSEnergy_HH
