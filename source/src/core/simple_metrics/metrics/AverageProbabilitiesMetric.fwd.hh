// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/AverageProbabilitiesMetric.fwd.hh
/// @brief A metric for averaging multiple PerResidueProbabilitiesMetrics
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

#ifndef INCLUDED_core_simple_metrics_metrics_AverageProbabilitiesMetric_fwd_hh
#define INCLUDED_core_simple_metrics_metrics_AverageProbabilitiesMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace core {
namespace simple_metrics {
namespace metrics {

class AverageProbabilitiesMetric;

using AverageProbabilitiesMetricOP = utility::pointer::shared_ptr< AverageProbabilitiesMetric >;
using AverageProbabilitiesMetricCOP = utility::pointer::shared_ptr< AverageProbabilitiesMetric const >;

} //metrics
} //simple_metrics
} //core

#endif //INCLUDED_core_simple_metrics_metrics_AverageProbabilitiesMetric_fwd_hh
