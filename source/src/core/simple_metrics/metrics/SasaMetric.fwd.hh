// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/metrics/SasaMetric.fwd.hh
/// @brief A Metric to cacluate overall sasa of a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_metrics_SasaMetric_fwd_hh
#define INCLUDED_core_simple_metrics_metrics_SasaMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace simple_metrics {
namespace metrics {

class SasaMetric;

typedef utility::pointer::shared_ptr< SasaMetric > SasaMetricOP;
typedef utility::pointer::shared_ptr< SasaMetric const > SasaMetricCOP;

} //core
} //simple_metrics
} //metrics

#endif //INCLUDED_core_metrics_SasaMetric_fwd_hh
