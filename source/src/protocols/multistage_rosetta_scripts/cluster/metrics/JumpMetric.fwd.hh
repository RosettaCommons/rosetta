// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/metrics/JumpMetric.fwd.hh
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_multistage_rosetta_scripts_cluster_metric_JumpMetric_fwd_hh
#define INCLUDED_protocols_multistage_rosetta_scripts_cluster_metric_JumpMetric_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {
namespace metrics {

class JumpMetric;
typedef utility::pointer::shared_ptr< JumpMetric > JumpMetricOP;
typedef utility::pointer::shared_ptr< JumpMetric const > JumpMetricCOP;

} // namespace metrics
} // namespace cluster
} // namespace multistage_rosetta_scripts
} // namespace protocols


#endif
