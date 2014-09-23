// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/NeighborsByDistanceCalculator.fwd.hh
/// @brief NeighborsByDistanceCalculator can determine all the neighbors of a residue within a certain distance.
/// @author Steven Lewis


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_NeighborsByDistanceCalculator_fwd_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_NeighborsByDistanceCalculator_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols{
namespace toolbox{
namespace pose_metric_calculators{

//Forwards and OP typedefs
class NeighborsByDistanceCalculator;
typedef utility::pointer::shared_ptr< NeighborsByDistanceCalculator > NeighborsByDistanceCalculatorOP;
typedef utility::pointer::shared_ptr< NeighborsByDistanceCalculator const > NeighborsByDistanceCalculatorCOP;

}//PoseMetricCalculators
}//toolbox
}//protocols

#endif //INCLUDED_protocols_toolbox_PoseMetricCalculators_NeighborsByDistanceCalculator_FWD_HH
