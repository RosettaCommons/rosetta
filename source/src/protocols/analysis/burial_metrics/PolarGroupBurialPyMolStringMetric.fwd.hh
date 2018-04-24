// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/burial_metrics/PolarGroupBurialPyMolStringMetric.fwd.hh
/// @brief Forward declarations for a string metric that generates a string of PyMol commands to colour
/// a structure's polar groups based on burial.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_analysis_burial_metrics_PolarGroupBurialPyMolStringMetric_fwd_hh
#define INCLUDED_protocols_analysis_burial_metrics_PolarGroupBurialPyMolStringMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace analysis {
namespace burial_metrics {

class PolarGroupBurialPyMolStringMetric;

typedef utility::pointer::shared_ptr< PolarGroupBurialPyMolStringMetric > PolarGroupBurialPyMolStringMetricOP;
typedef utility::pointer::shared_ptr< PolarGroupBurialPyMolStringMetric const > PolarGroupBurialPyMolStringMetricCOP;

} //protocols
} //analysis
} //burial_metrics

#endif //INCLUDED_protocols_analysis_burial_metrics_PolarGroupBurialPyMolStringMetric_fwd_hh
