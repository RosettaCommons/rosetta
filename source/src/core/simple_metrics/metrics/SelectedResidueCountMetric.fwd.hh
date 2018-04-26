// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SelectedResidueCountMetric.fwd.hh
/// @brief A SimpleMetric that counts the number of residues in a residue selection.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_core_simple_metrics_metrics_SelectedResidueCountMetric_fwd_hh
#define INCLUDED_core_simple_metrics_metrics_SelectedResidueCountMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace simple_metrics {
namespace metrics {

class SelectedResidueCountMetric;

typedef utility::pointer::shared_ptr< SelectedResidueCountMetric > SelectedResidueCountMetricOP;
typedef utility::pointer::shared_ptr< SelectedResidueCountMetric const > SelectedResidueCountMetricCOP;

} //protocols
} //analysis
} //simple_metrics

#endif //INCLUDED_core_simple_metrics_metrics_SelectedResidueCountMetric_fwd_hh
