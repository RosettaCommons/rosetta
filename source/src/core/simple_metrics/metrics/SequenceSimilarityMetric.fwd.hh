// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SequenceSimilarityMetric.fwd.hh
/// @brief compares the sequences of the native structure (using native flag) to the sequence of a given pose.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_simple_metrics_metrics_SequenceSimilarityMetric_fwd_hh
#define INCLUDED_core_simple_metrics_metrics_SequenceSimilarityMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace core {
namespace simple_metrics {
namespace metrics {

class SequenceSimilarityMetric;

typedef utility::pointer::shared_ptr< SequenceSimilarityMetric > SequenceSimilarityMetricOP;
typedef utility::pointer::shared_ptr< SequenceSimilarityMetric const > SequenceSimilarityMetricCOP;

} //core
} //simple_metrics
} //metrics

#endif //INCLUDED_core_metrics_SequenceSimilarityMetric_fwd_hh
