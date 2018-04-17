// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/CompositeStringMetric.hh
///
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_test_classes_fwd_hh
#define INCLUDED_core_simple_metrics_test_classes_fwd_hh

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

namespace core {
namespace simple_metrics {

class TestStringMetric;
typedef utility::pointer::shared_ptr< TestStringMetric > TestStringMetricOP;
typedef utility::pointer::shared_ptr< TestStringMetric const > TestStringMetricCOP;

class TestRealMetric;
typedef utility::pointer::shared_ptr< TestRealMetric > TestRealMetricOP;
typedef utility::pointer::shared_ptr< TestRealMetric const > TestRealMetricCOP;

class TestCompositeStringMetric;
typedef utility::pointer::shared_ptr< TestCompositeStringMetric > TestCompositeStringMetricOP;
typedef utility::pointer::shared_ptr< TestCompositeStringMetric const > TestCompositeStringMetricCOP;


class TestCompositeRealMetric;
typedef utility::pointer::shared_ptr< TestCompositeRealMetric > TestCompositeRealMetricOP;
typedef utility::pointer::shared_ptr< TestCompositeRealMetric const > TestCompositeRealMetricCOP;

}
}

#endif
