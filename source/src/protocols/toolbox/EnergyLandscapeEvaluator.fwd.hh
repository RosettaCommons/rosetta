// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/EnergyLandscapeEvaluator.fwd.hh
/// @brief Evaluates a set of score/rms points
/// @author Tom Linsky (tlinsky@gmail.com)


#ifndef INCLUDED_protocols_toolbox_EnergyLandscapeEvaluator_fwd_hh
#define INCLUDED_protocols_toolbox_EnergyLandscapeEvaluator_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace toolbox {

class ScoreRmsPoint;
typedef utility::pointer::shared_ptr< ScoreRmsPoint > ScoreRmsPointOP;
typedef utility::pointer::shared_ptr< ScoreRmsPoint const > ScoreRmsPointCOP;

class ScoreRmsPoints;
typedef utility::pointer::shared_ptr< ScoreRmsPoints > ScoreRmsPointsOP;
typedef utility::pointer::shared_ptr< ScoreRmsPoints const > ScoreRmsPointsCOP;

class EnergyLandscapeEvaluator;
typedef utility::pointer::shared_ptr< EnergyLandscapeEvaluator > EnergyLandscapeEvaluatorOP;
typedef utility::pointer::shared_ptr< EnergyLandscapeEvaluator const > EnergyLandscapeEvaluatorCOP;

} //protocols
} //toolbox


#endif //INCLUDED_protocols_toolbox_EnergyLandscapeEvaluator_fwd_hh

