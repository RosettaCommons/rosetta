// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/InterGroupNeighborsCalculator.fwd.hh
/// @brief InterGroupNeighborsCalculator can determine all the neighbors of a residue within a certain distance.
/// @author Steven Lewis


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_InterGroupNeighborsCalculator_fwd_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_InterGroupNeighborsCalculator_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace pose_metric_calculators {

//Forwards and OP typedefs
class InterGroupNeighborsCalculator;
typedef utility::pointer::shared_ptr< InterGroupNeighborsCalculator > InterGroupNeighborsCalculatorOP;
typedef utility::pointer::shared_ptr< InterGroupNeighborsCalculator const > InterGroupNeighborsCalculatorCOP;

}//PoseMetricCalculators
}//protocols

#endif //INCLUDED_protocols_toolbox_PoseMetricCalculators_InterGroupNeighborsCalculator_FWD_HH
