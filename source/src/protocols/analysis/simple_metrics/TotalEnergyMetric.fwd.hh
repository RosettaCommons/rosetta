// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/TotalEnergyMetric.fwd.hh
/// @brief A metric to report the total energy of the system or the delta total energy between another input pose or the set native.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_analysis_simple_metrics_TotalEnergyMetric_fwd_hh
#define INCLUDED_protocols_analysis_simple_metrics_TotalEnergyMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace analysis {
namespace simple_metrics {

class TotalEnergyMetric;

typedef utility::pointer::shared_ptr< TotalEnergyMetric > TotalEnergyMetricOP;
typedef utility::pointer::shared_ptr< TotalEnergyMetric const > TotalEnergyMetricCOP;

} //core
} //simple_metrics
} //metrics

#endif //INCLUDED_protocols_analysis_simple_metrics_TotalEnergyMetric_fwd_hh
