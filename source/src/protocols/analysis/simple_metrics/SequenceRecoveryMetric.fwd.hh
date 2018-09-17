// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/simple_metrics/SequenceRecoveryMetric.fwd.hh
/// @brief Calculate sequence recovery statistics on a protein, relative to a reference.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_analysis_simple_metrics_SequenceRecoveryMetric_fwd_hh
#define INCLUDED_protocols_analysis_simple_metrics_SequenceRecoveryMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace analysis {
namespace simple_metrics {

class SequenceRecoveryMetric;

typedef utility::pointer::shared_ptr< SequenceRecoveryMetric > SequenceRecoveryMetricOP;
typedef utility::pointer::shared_ptr< SequenceRecoveryMetric const > SequenceRecoveryMetricCOP;

} //protocols
} //analysis
} //simple_metrics

#endif //INCLUDED_protocols_analysis_simple_metrics_SequenceRecoveryMetric_fwd_hh
