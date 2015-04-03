// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh
/// @brief  Context-Dependent, Long-Range, Two-Body Energy class declaration
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_methods_ContextDependentLRTwoBodyEnergy_hh
#define INCLUDED_core_scoring_methods_ContextDependentLRTwoBodyEnergy_hh

// Unit headers
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

class ContextDependentLRTwoBodyEnergy : public LongRangeTwoBodyEnergy {
public:
	typedef LongRangeTwoBodyEnergy parent;

public:
	/// @brief Constructor to inform the EnergyMethod ancestor which
	/// ScoreTypes this EnergyMethod is responsible for computing
	ContextDependentLRTwoBodyEnergy( EnergyMethodCreatorOP );

	virtual ~ContextDependentLRTwoBodyEnergy();


	EnergyMethodType
	method_type() const
	{
		return cd_lr_2b;
	}

};


} // methods
} // scoring
} // core


#endif

