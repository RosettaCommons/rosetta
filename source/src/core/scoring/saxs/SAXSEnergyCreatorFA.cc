// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/saxs/SAXSEnergy.cc
/// @brief  "Energy" based on a similarity of theoretical SAXS spectrum computed for a pose and the experimental data
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <core/scoring/saxs/SAXSEnergyFA.hh>
#include <core/scoring/saxs/SAXSEnergyCreatorFA.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace saxs {


ScoreTypes SAXSEnergyCreatorFA::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( saxs_fa_score );
	return sts;
}

methods::EnergyMethodOP SAXSEnergyCreatorFA::create_energy_method( methods::EnergyMethodOptions const &) const {
	return methods::EnergyMethodOP( new SAXSEnergyFA() );
}
core::Size
SAXSEnergyCreatorFA::version() const
{
	return 1; // Initial versioning
}

} // saxs
} // scoring
} // core
