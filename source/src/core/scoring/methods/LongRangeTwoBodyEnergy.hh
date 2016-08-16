// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/LongRangeTwoBodyEnergy.hh
/// @brief  Long range two body energy method class declaration
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_methods_LongRangeTwoBodyEnergy_hh
#define INCLUDED_core_scoring_methods_LongRangeTwoBodyEnergy_hh

// Unit Headers
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>

// Package Headers
#include <core/scoring/methods/TwoBodyEnergy.hh>
#include <core/scoring/methods/Methods.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

class LongRangeTwoBodyEnergy : public TwoBodyEnergy
{
public:
	typedef TwoBodyEnergy parent;

public:
	/// @brief Constructor with an EnergyMethodCreator to inform the grandparent EnergyMethod
	/// class which ScoreTypes this EnergyMethod computes.
	LongRangeTwoBodyEnergy( EnergyMethodCreatorOP creator );

	virtual ~LongRangeTwoBodyEnergy();

	virtual LongRangeEnergyType long_range_type() const = 0;

	virtual
	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const = 0;

};

}
}
}

#endif
