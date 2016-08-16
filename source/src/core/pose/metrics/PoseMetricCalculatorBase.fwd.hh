// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author John Karanicolas


#ifndef INCLUDED_core_pose_metrics_PoseMetricCalculatorBase_fwd_hh
#define INCLUDED_core_pose_metrics_PoseMetricCalculatorBase_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pose {
namespace metrics {

class PoseMetricCalculator;
typedef utility::pointer::shared_ptr< PoseMetricCalculator > PoseMetricCalculatorOP;
typedef utility::pointer::shared_ptr< PoseMetricCalculator const > PoseMetricCalculatorCOP;

} // namespace metrics
} // namespace pose
} // namespace core


#endif
