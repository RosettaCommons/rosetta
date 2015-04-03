// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_methods_ContextDependentOneBodyEnergy_hh
#define INCLUDED_core_scoring_methods_ContextDependentOneBodyEnergy_hh

// Unit headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/OneBodyEnergy.hh>

#include <core/scoring/EnergyMap.fwd.hh>


// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class ContextDependentOneBodyEnergy : public OneBodyEnergy {
public:
	typedef OneBodyEnergy parent;

public:
	/// @brief Constructor with an EnergyMethodCreator to inform the EnergyMethod
	/// parent which ScoreTypes this EnergyMethod is responsible for computing.
	ContextDependentOneBodyEnergy( EnergyMethodCreatorOP );

	/// @brief Returns the cd_1b element of the EnergyMethodType enumeration; this method
	/// should NOT be overridden by derived classes.
	virtual
	EnergyMethodType
	method_type() const;


	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const = 0;

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
