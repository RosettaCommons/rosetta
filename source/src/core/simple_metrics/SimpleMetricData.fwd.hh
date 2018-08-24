// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/SimpleMetricData.fwd.hh
/// @brief
/// @author Jared Adolf-Bryfogle


#ifndef INCLUDED_basic_datacache_SimpleMetricData_fwd_hh
#define INCLUDED_basic_datacache_SimpleMetricData_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/down_cast.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>


namespace core {
namespace simple_metrics {


class SimpleMetricData;
typedef utility::pointer::shared_ptr< SimpleMetricData > SimpleMetricDataOP;
typedef utility::pointer::shared_ptr< SimpleMetricData const > SimpleMetricDataCOP;
typedef utility::pointer::weak_ptr< SimpleMetricData > SimpleMetricDataAP;
typedef utility::pointer::weak_ptr< SimpleMetricData const > SimpleMetricDataCAP;


} // namespace simple_metrics
} // namespace core


#endif /* INCLUDED_core_simple_metrics_SimpleMetricData_FWD_HH */
