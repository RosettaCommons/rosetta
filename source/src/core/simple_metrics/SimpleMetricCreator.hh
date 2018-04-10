// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/SimpleMetricCreator.hh
/// @brief  Class for instantiating a particular SimpleMetric
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_SimpleMetricCreator_HH
#define INCLUDED_core_simple_metrics_SimpleMetricCreator_HH

// Package headers
#include <core/simple_metrics/SimpleMetric.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace core {
namespace simple_metrics {

class SimpleMetricCreator : public utility::pointer::ReferenceCount {
public:
	/// @brief Instantiate a particular SimpleMetric
	virtual SimpleMetricOP
	create_simple_metric() const = 0;

	/// @brief Return a string that will be used to instantiate the particular SimpleMetric
	virtual std::string
	keyname() const = 0;

	virtual void
	provide_xml_schema( utility::tag::XMLSchemaDefinition &) const = 0;
};


} //namespace simple_metrics
} //namespace core


#endif
