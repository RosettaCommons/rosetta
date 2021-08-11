// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/simple_metrics/BlockwisePoseCompMotifMetric.fwd.hh
/// @brief composite metric returning motif scores that pass per block
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_pose_sewing_simple_metrics_BlockwisePoseCompMotifMetric_fwd_hh
#define INCLUDED_protocols_pose_sewing_simple_metrics_BlockwisePoseCompMotifMetric_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace pose_sewing {
namespace simple_metrics {

class BlockwisePoseCompMotifMetric;

using BlockwisePoseCompMotifMetricOP = utility::pointer::shared_ptr< BlockwisePoseCompMotifMetric >;
using BlockwisePoseCompMotifMetricCOP = utility::pointer::shared_ptr< BlockwisePoseCompMotifMetric const >;

} //simple_metrics
} //pose_sewing
} //protocols

#endif //INCLUDED_protocols_pose_sewing_simple_metrics_BlockwisePoseCompMotifMetric_fwd_hh
