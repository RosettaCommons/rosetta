// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/analysis/simple_metrics/SimpleMetricRegistrator.hh
/// @brief  Declaration of the template class for registrating SimpleMetricCreators with
///         the SimpleMetricFactory
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_SimpleMetricRegistrator_hh
#define INCLUDED_core_simple_metrics_SimpleMetricRegistrator_hh

// Package headers
#include <core/simple_metrics/SimpleMetricFactory.hh>
#include <utility/factory/WidgetRegistrator.hh>

namespace core {
namespace simple_metrics {


/// @brief This templated class will register an instance of an
/// SimpleMetricCreator (class T) with the SimpleMetricFactory.  It will ensure
/// that no SimpleMetric creator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class SimpleMetricRegistrator : public utility::factory::WidgetRegistrator< SimpleMetricFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< SimpleMetricFactory, T > parent;
public:
	SimpleMetricRegistrator() : parent() {}
};

}
}

#endif
