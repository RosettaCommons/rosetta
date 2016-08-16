// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ContextIndependentTwoBodyEnergy.hh
/// @brief  Short ranged, context-independent, two-body energy class declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_methods_ContextIndependentTwoBodyEnergy_hh
#define INCLUDED_core_scoring_methods_ContextIndependentTwoBodyEnergy_hh

// Unit headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

class ContextIndependentTwoBodyEnergy : public ShortRangeTwoBodyEnergy {
public:
	typedef ShortRangeTwoBodyEnergy parent;

public:
	/// @brief Constructor with an EnergyMethodCreator to inform the
	/// ancestor EnergyMethod class which ScoreTypes this EnergyMethod
	/// is responsible for computing.
	ContextIndependentTwoBodyEnergy( EnergyMethodCreatorOP );

	virtual
	~ContextIndependentTwoBodyEnergy();

	EnergyMethodType
	method_type() const;

};

} // methods
} // scoring
} // core

#endif
