// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/SimpleMetricFactory.hh
/// @brief  Class for instantiating arbitrary SimpleMetrics from a string --> SimpleMetricCreator map
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_SimpleMetricFactory_HH
#define INCLUDED_core_simple_metrics_SimpleMetricFactory_HH

// Unit headers
#include <utility/SingletonBase.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
// Package headers
#include <core/simple_metrics/SimpleMetric.fwd.hh>
#include <core/simple_metrics/SimpleMetricCreator.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>

// C++ headers
#include <map>
#include <string>

namespace core {
namespace simple_metrics {

class SimpleMetricFactory : public utility::SingletonBase< SimpleMetricFactory > {
private:
	typedef std::map< std::string, SimpleMetricCreatorOP > CreatorMap;

public:
	SimpleMetricFactory();

	void factory_register( SimpleMetricCreatorOP creator );

	bool has_type( std::string const & simple_metric_name ) const;

	SimpleMetricOP new_simple_metric(
		std::string const & constraint_generator_name,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) const;

	/// @brief Get the XML schema for a given residue selector.
	/// @details Throws an error if the residue selector is unknown to Rosetta.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void provide_xml_schema( std::string const &selector_name, utility::tag::XMLSchemaDefinition & xsd ) const;

	void define_simple_metric_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

	//All of these may move to a new file for schema generation

	static std::string simple_metric_xml_schema_group_name();

private:
	CreatorMap creator_map_;
};

} //namespace simple_metrics
} //namespace core


#endif
