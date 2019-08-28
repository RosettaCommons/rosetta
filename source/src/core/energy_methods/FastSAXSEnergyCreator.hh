// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/saxs/FastSAXSEnergyCreator.hh
/// @brief  Energy Creator for FastSAX scoring of Stovgaard et al (BMC Bioinf. 2010)
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_saxs_FastSAXSEnergyCreator_hh
#define INCLUDED_core_scoring_saxs_FastSAXSEnergyCreator_hh

#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/types.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace saxs {

class FastSAXSEnergyCreator : public methods::EnergyMethodCreator {
public:

	FastSAXSEnergyCreator() {}

	/// @brief Instantiate a new SAXSEnergy
	virtual methods::EnergyMethodOP create_energy_method(methods::EnergyMethodOptions const &) const;

	/// @brief Return the set of score types claimed by the EnergyMethod
	/// this EnergyMethodCreator creates in its create_energy_method() function
	virtual ScoreTypes score_types_for_method() const;
};

}
}
}

#endif
