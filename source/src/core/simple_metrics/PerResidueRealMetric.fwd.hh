// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/PerResidueRealMetric.hh
///
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_PerResidueRealMetric_fwd_hh
#define INCLUDED_core_simple_metrics_PerResidueRealMetric_fwd_hh

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

namespace core {
namespace simple_metrics {

class PerResidueRealMetric;
typedef utility::pointer::shared_ptr< PerResidueRealMetric > PerResidueRealMetricOP;
typedef utility::pointer::shared_ptr< PerResidueRealMetric const > PerResidueRealMetricCOP;

}
}

#endif
