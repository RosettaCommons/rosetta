// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/methods/SplitUnfoldedTwoBodyEnergyCreator.hh
/// @brief  Energy creator for the split unfolded two body energy method
/// @author Riley Simmons-Edler (rse231@nyu.edu)

#ifndef INCLUDED_core_scoring_methods_SplitUnfoldedTwoBodyEnergyCreator_hh
#define INCLUDED_core_scoring_methods_SplitUnfoldedTwoBodyEnergyCreator_hh

#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>

#include <utility/vector1.hh>

#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace scoring {
namespace methods {

class SplitUnfoldedTwoBodyEnergyCreator : public EnergyMethodCreator
{
public:
	virtual methods::EnergyMethodOP create_energy_method(const methods::EnergyMethodOptions &) const;

	/// @brief Return the set of score types claimed by the EnergyMethod
	/// this EnergyMethodCreator creates in its create_energy_method() function
	virtual ScoreTypes score_types_for_method() const;
};

}
}
}

#endif
