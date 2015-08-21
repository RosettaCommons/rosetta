// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/vall_lookback/VallLookbackEnergy.hh
/// @brief  lookback to native scores.
/// @author TJ Brunette (brunette@gmail.com)

/// This energy function requires A. the VALL lookback database gets added to the potential
/// B. For efficiency when locations are sampled they need to be marked changed in the VallLookback datacache object

#include <core/scoring/methods/vall_lookback/VallLookbackEnergy.hh>
#include <core/scoring/methods/vall_lookback/VallLookbackEnergyCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/conformation/Residue.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/indexed_structure_store.OptionKeys.gen.hh>

namespace core {
namespace scoring {
namespace methods {

methods::EnergyMethodOP
VallLookbackEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP(new VallLookbackEnergy);
}

ScoreTypes
VallLookbackEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( vall_lb );
	return sts;
}

/// c-tor
VallLookbackEnergy::VallLookbackEnergy() :
	parent( methods::EnergyMethodCreatorOP(new VallLookbackEnergyCreator )),
	potential_( ScoringManager::get_instance()->get_vallLookbackPotential())
{
}

/// clone
EnergyMethodOP
VallLookbackEnergy::clone() const
{
	return EnergyMethodOP( new VallLookbackEnergy);
}

void
VallLookbackEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	using namespace basic::options;
	using namespace OptionKeys::indexed_structure_store;
	using namespace conformation;
	Real maxRmsd = potential_.lookback(pose);
	Real thresholdDistance = option[fragment_threshold_distance]();
	if ( maxRmsd <= thresholdDistance ) {
		totals[vall_lb] = 0;
	} else {
		totals[vall_lb] = maxRmsd -thresholdDistance;
	}

} // finalize_total_energy

core::Size
VallLookbackEnergy::version() const {
	return 1; // Initial version
}


} // methods
} // scoring
} // core
