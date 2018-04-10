// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/metrics/Metric.fwd.hh
/// @brief The base class for Metrics in the Metric/Filter/Reporter system
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_SimpleMetric_fwd_hh
#define INCLUDED_core_simple_metrics_SimpleMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace simple_metrics {

class SimpleMetric;

typedef utility::pointer::shared_ptr< SimpleMetric > SimpleMetricOP;
typedef utility::pointer::shared_ptr< SimpleMetric const > SimpleMetricCOP;

} //core
} //metrics

#endif //INCLUDED_core_metrics_SimpleMetric_fwd_hh
