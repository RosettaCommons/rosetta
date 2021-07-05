// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreFunction.hh
/// @brief  Score function class
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_methods_ContextDependentTwoBodyEnergy_hh
#define INCLUDED_core_scoring_methods_ContextDependentTwoBodyEnergy_hh

// Unit headers
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>

// Project headers
//#include <core/pack/rotamer_set/RotamerSet.fwd.hh>

// ObjexxFCL headers



namespace core {
namespace scoring {
namespace methods {


class ContextDependentTwoBodyEnergy : public ShortRangeTwoBodyEnergy {
public:
	typedef ShortRangeTwoBodyEnergy parent;

public:
	ContextDependentTwoBodyEnergy( EnergyMethodCreatorOP creator );

	~ContextDependentTwoBodyEnergy() override;


	EnergyMethodType
	method_type() const override;

};

} // methods
} // scoring
} // core


#endif
