// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/UnfoldedStateEnergy.hh
/// @brief  Unfolded state energy declaration header file
/// @author Ron Jacak (ronj@email.unc.edu)


#ifndef INCLUDED_core_energy_methods_UnfoldedStateEnergy_hh
#define INCLUDED_core_energy_methods_UnfoldedStateEnergy_hh

// Unit headers
#include <core/energy_methods/UnfoldedStateEnergy.fwd.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>

#include <core/scoring/UnfoldedStatePotential.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



class UnfoldedStateEnergy : public core::scoring::methods::ContextIndependentOneBodyEnergy {

public:
	typedef core::scoring::methods::ContextIndependentOneBodyEnergy parent;

	UnfoldedStateEnergy( std::string const & type );
	UnfoldedStateEnergy( std::string const & type, const core::scoring::EnergyMap & emap_in );
	~UnfoldedStateEnergy() override;

	core::scoring::methods::EnergyMethodOP clone() const override;

	void
	residue_energy( conformation::Residue const & rsd, pose::Pose const & pose, core::scoring::EnergyMap & emap ) const override;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override { return false; }

	void indicate_required_context_graphs( utility::vector1< bool > & ) const override;

private:

	std::string type_;
	core::scoring::UnfoldedStatePotential const & unf_state_potential_;
	core::scoring::EnergyMap score_type_weights_;
	core::Size version() const override;

};

} // scoring
} // core


#endif // INCLUDED_core_energy_methods_UnfoldedStateEnergy_HH
