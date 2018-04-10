// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/TimingProfileMetric.fwd.hh
/// @brief Calculate the time difference between construction and apply/calculate.  Useful to time protocols in RosettaScripts or through mover containers.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_metrics_TimingProfileMetric_fwd_hh
#define INCLUDED_core_simple_metrics_metrics_TimingProfileMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace simple_metrics {
namespace metrics {

class TimingProfileMetric;

typedef utility::pointer::shared_ptr< TimingProfileMetric > TimingProfileMetricOP;
typedef utility::pointer::shared_ptr< TimingProfileMetric const > TimingProfileMetricCOP;

} //core
} //simple_metrics
} //metrics

#endif //INCLUDED_core_simple_metrics_metrics_TimingProfileMetric_fwd_hh
