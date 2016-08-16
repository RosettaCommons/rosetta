// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh
/// @brief  Context-Independent, Long-Range, Two-Body Energy class declaration
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_methods_ContextIndependentLRTwoBodyEnergy_hh
#define INCLUDED_core_scoring_methods_ContextIndependentLRTwoBodyEnergy_hh

// Unit headers
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

class ContextIndependentLRTwoBodyEnergy : public LongRangeTwoBodyEnergy {
public:
	typedef LongRangeTwoBodyEnergy parent;

public:
	ContextIndependentLRTwoBodyEnergy( EnergyMethodCreatorOP );

	virtual ~ContextIndependentLRTwoBodyEnergy();


	EnergyMethodType
	method_type() const;

};


} // methods
} // scoring
} // core


#endif
