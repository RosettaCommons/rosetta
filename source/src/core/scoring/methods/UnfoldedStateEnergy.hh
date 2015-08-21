// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/UnfoldedStateEnergy.hh
/// @brief  Unfolded state energy declaration header file
/// @author Ron Jacak (ronj@email.unc.edu)


#ifndef INCLUDED_core_scoring_methods_UnfoldedStateEnergy_hh
#define INCLUDED_core_scoring_methods_UnfoldedStateEnergy_hh

// Unit headers
#include <core/scoring/methods/UnfoldedStateEnergy.fwd.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/UnfoldedStatePotential.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class UnfoldedStateEnergy : public ContextIndependentOneBodyEnergy {

public:
	typedef ContextIndependentOneBodyEnergy parent;

	UnfoldedStateEnergy( std::string const & type );
	UnfoldedStateEnergy( std::string const & type, const EnergyMap & emap_in );
	~UnfoldedStateEnergy();

	virtual
	EnergyMethodOP clone() const;

	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap ) const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;

private:

	std::string type_;
	UnfoldedStatePotential const & unf_state_potential_;
	EnergyMap score_type_weights_;
	virtual
	core::Size version() const;

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_methods_UnfoldedStateEnergy_HH
