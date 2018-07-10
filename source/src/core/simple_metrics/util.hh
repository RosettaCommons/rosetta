// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/util.hh
/// @brief Util files for SimpleMetrics.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_util_hh
#define INCLUDED_core_simple_metrics_util_hh

#include <core/simple_metrics/SimpleMetric.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

//C++ headers
#include <map>

namespace core {
namespace simple_metrics {

std::string
complex_type_name_for_simple_metric( std::string const & metric_name );

///@brief Generate the ComplexTypeGenerator from the SimpleMetric base class.
///  Add any additional schema options from sub-derived classes
void
xsd_simple_metric_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & simple_metric_name,
	std::string const & description,
	utility::tag::AttributeList const & attributes);

///@brief Generate the ComplexTypeGenerator from the SimpleMetric and PerResidueRealMetric base classes.
///  Add any additional schema options from sub-derived classes
void
xsd_per_residue_real_metric_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & simple_metric_name,
	std::string const & description,
	utility::tag::AttributeList const & attributes);

///@brief Generate the ComplexTypeGenerator from the SimpleMetric and PerResidueStringMetric base classes.
///  Add any additional schema options from sub-derived classes
void
xsd_per_residue_string_metric_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & simple_metric_name,
	std::string const & description,
	utility::tag::AttributeList const & attributes);

utility::vector1< SimpleMetricCOP >
get_metrics_from_datamap_and_subtags(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	std::string tag_name="metrics");

SimpleMetricCOP
get_metric_from_datamap_and_subtags(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	std::string tag_name="metric");

///@brief Add the schema for per-residue simple metrics.  This is to reduce code dup.
void
add_per_residue_simple_metric_schema( utility::tag::XMLSchemaComplexTypeGeneratorOP complex_schema );

} //core
} //simple_metrics


#endif //core/simple_metrics_util_hh

