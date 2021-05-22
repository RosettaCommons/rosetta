// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/simple_metrics/ConstraintsMetric.fwd.hh
/// @brief A simple metric that writes out all the constraints in a pose or sub-region of a pose,
/// in a format matching the (non-enzdes) constraints file format.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_analysis_simple_metrics_ConstraintsMetric_fwd_hh
#define INCLUDED_protocols_analysis_simple_metrics_ConstraintsMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace analysis {
namespace simple_metrics {

class ConstraintsMetric;

using ConstraintsMetricOP = utility::pointer::shared_ptr< ConstraintsMetric >;
using ConstraintsMetricCOP = utility::pointer::shared_ptr< ConstraintsMetric const >;

} //simple_metrics
} //analysis
} //protocols

#endif //INCLUDED_protocols_analysis_simple_metrics_ConstraintsMetric_fwd_hh
